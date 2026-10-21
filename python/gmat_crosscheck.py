"""Cross-check against NASA GMAT.

1. GMAT's differential corrector solves the same Hohmann transfers that odst.py sizes
   analytically, and the GTO -> GEO apogee burn of the apse-to-apse model.
2. GMAT flies every elliptical-mission burn plan (apse-to-apse and Lambert) from its
   start state, and the orbit it arrives on is compared with the target orbit.

Needs a local GMAT (https://sourceforge.net/projects/gmat/). Not run in CI.
Run: python python/gmat_crosscheck.py [path/to/GmatConsole.exe]
     (default: $GMAT_CONSOLE, else D:/GMAT/bin/GmatConsole.exe)
"""
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from string import Template

from odst import MU, RE, Mission, evaluate, hohmann, target_mismatch, vcirc, vvis

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
TEMPLATE = Template((HERE.parent / "gmat" / "hohmann_plane_change.script").read_text())
APOGEE = Template((HERE.parent / "gmat" / "apogee_burn.script").read_text())
POINT_MASS = "FM.PointMasses = {Earth};"
J2 = ("FM.PrimaryBodies = {Earth};\nFM.GravityField.Earth.Degree = 2;\n"
      "FM.GravityField.Earth.Order = 0;\nFM.GravityField.Earth.PotentialFile = 'JGM2.cof';")
INC1 = 30.0   # start inclination; target is INC1 - di (see the scripts' headers)
ELLIPTICAL = [Mission(alt1=250, apo1=35786, alt2=35786, di_deg=28.5),                   # GTO -> GEO
              Mission(alt1=400, alt2=600, apo2=39700, dw_deg=270, di_deg=0),            # LEO -> Molniya
              Mission(alt1=500, apo1=5000, alt2=800, apo2=12000, dw_deg=60, di_deg=0),  # apse lines 60 deg apart
              Mission(alt1=500, apo1=5000, alt2=800, apo2=12000, dw_deg=60, di_deg=10)]


def gmat(console, make_script):
    """Run the script make_script(report_path) builds; return the reported numbers or None."""
    with tempfile.TemporaryDirectory() as tmp:
        report, script = Path(tmp) / "report.txt", Path(tmp) / "case.script"
        script.write_text(make_script(report.as_posix()))
        subprocess.run([console, str(script)], cwd=Path(console).parent,
                       capture_output=True, text=True, timeout=600)
        return [float(x) for x in report.read_text().split()] if report.exists() else None


def solve_hohmann(console, alt1, alt2, di_deg, gravity):
    """The normal component's sign depends on the node; try both, keep the one that converged."""
    r1, r2, di = RE + alt1, RE + alt2, math.radians(di_deg)
    vt1, vt2 = vvis(r1, (r1 + r2) / 2), vvis(r2, (r1 + r2) / 2)
    for sign in (1, -1):
        v = gmat(console, lambda rep: TEMPLATE.substitute(
            r1=f"{r1:.6f}", r2=f"{r2:.6f}", inc1=INC1, inc2=INC1 - di_deg, gravity=gravity,
            dv1=f"{vt1 - vcirc(r1):.6f}", dv2v=f"{vcirc(r2) * math.cos(di) - vt2:.6f}",
            dv2n=f"{sign * vcirc(r2) * math.sin(di):.6f}", report=rep))
        if v and abs(v[9] - (INC1 - di_deg)) < 1e-6 and abs(v[7] - r2) < 1e-4:
            return math.hypot(*v[0:3]) + math.hypot(*v[3:6]), v[6]
    return None


def solve_apogee(console, m):
    """GMAT targets the single apoapsis burn from orbit (m.r1, m.ra1) onto a circle of radius m.ra1."""
    a, e, r2, di = (m.r1 + m.ra1) / 2, (m.ra1 - m.r1) / (m.ra1 + m.r1), m.ra1, math.radians(m.di_deg)
    va = vvis(r2, a)
    for sign in (1, -1):
        v = gmat(console, lambda rep: APOGEE.substitute(
            sma=f"{a:.6f}", ecc=f"{e:.12f}", inc1=INC1, inc2=INC1 - m.di_deg, r2=f"{r2:.6f}",
            dvv=f"{vcirc(r2) * math.cos(di) - va:.6f}", dvn=f"{sign * vcirc(r2) * math.sin(di):.6f}", report=rep))
        if v and abs(v[5] - (INC1 - m.di_deg)) < 1e-6 and abs(v[3] - r2) < 1e-4:
            return math.hypot(*v[0:3])
    return None


def fly_script(y0, plan, report):
    """GMAT script that starts from Cartesian y0 and applies inertial burns at plan times."""
    s = ["Create Spacecraft Sat;", "Sat.DateFormat = UTCGregorian;", "Sat.Epoch = '01 Jan 2026 12:00:00.000';",
         "Sat.CoordinateSystem = EarthMJ2000Eq;", "Sat.DisplayStateType = Cartesian;"]
    s += [f"Sat.{k} = {x:.17g};" for k, x in zip(("X", "Y", "Z", "VX", "VY", "VZ"), y0)]
    s += ["Create ForceModel FM;", "FM.CentralBody = Earth;", POINT_MASS, "Create Propagator Prop;", "Prop.FM = FM;",
          "Prop.Type = RungeKutta89;", "Prop.InitialStepSize = 60;", "Prop.Accuracy = 1e-13;", "Prop.MinStep = 0;",
          "Prop.MaxStep = 600;"]
    for k, (_, dv) in enumerate(plan):
        s += [f"Create ImpulsiveBurn B{k};", f"B{k}.CoordinateSystem = EarthMJ2000Eq;"]
        s += [f"B{k}.Element{j + 1} = {dv[j]:.17g};" for j in range(3)]
    s += ["Create ReportFile Rpt;", f"Rpt.Filename = '{report}';", "Rpt.WriteHeaders = false;", "Rpt.Precision = 17;",
          "BeginMissionSequence;"]
    t = 0.0
    for k, (tb, _) in enumerate(plan):
        if tb > t:
            s.append(f"Propagate Prop(Sat) {{Sat.ElapsedSecs = {tb - t:.17g}}};")
        s.append(f"Maneuver B{k}(Sat);")
        t = tb
    s += ["Propagate Prop(Sat) {Sat.ElapsedSecs = 3600};",
          "Report Rpt " + " ".join(f"Sat.EarthMJ2000Eq.{k}" for k in ("X", "Y", "Z", "VX", "VY", "VZ")) + ";"]
    return "\n".join(s) + "\n"


def main():
    console = (sys.argv[1] if len(sys.argv) > 1 else
               os.environ.get("GMAT_CONSOLE", "D:/GMAT/bin/GmatConsole.exe"))
    fails = 0

    def report(ok, msg):
        nonlocal fails
        fails += not ok
        print(f"  [{'PASS' if ok else 'FAIL'}] {msg}")

    print(f"GMAT: {console}\n\n=== 1. GMAT targets the transfers itself ===")
    print(f"{'Case':40s}{'ODST dv':>10s}{'GMAT dv':>10s}{'diff m/s':>10s}{'ODST TOF':>10s}{'GMAT TOF':>10s}{'diff s':>9s}")
    for name, alt1, alt2, di, grav in [("400 km -> GEO, di 28.5 deg, two-body", 400, 35786, 28.5, POINT_MASS),
                                       ("400 km -> GEO, coplanar, two-body", 400, 35786, 0.0, POINT_MASS),
                                       ("400 km -> GEO, di 28.5 deg, J2 (JGM-2)", 400, 35786, 28.5, J2)]:
        dv, tof = hohmann(RE + alt1, RE + alt2, math.radians(di))[:2]
        g = solve_hohmann(console, alt1, alt2, di, grav)
        if g is None:
            report(False, f"{name}: GMAT did not converge")
            continue
        ddv, dt = 1e3 * (g[0] - dv), g[1] - tof
        print(f"{name:40s}{dv:10.5f}{g[0]:10.5f}{ddv:10.4f}{tof / 3600:9.4f}h{g[1] / 3600:9.4f}h{dt:9.3f}")
        if grav == POINT_MASS:        # same physics: must agree to 1 mm/s and 0.01 s
            report(abs(ddv) < 1e-3 and abs(dt) < 0.01, "two-body agreement within 1 mm/s and 0.01 s")
    gto = ELLIPTICAL[0]
    apse = next(o for o in evaluate(gto) if o.key == "apse")
    g = solve_apogee(console, gto)
    if g is None:
        report(False, "GTO -> GEO apogee burn: GMAT did not converge")
    else:
        print(f"{'GTO 250 x 35786 km -> GEO, di 28.5 deg':40s}{apse.dv:10.5f}{g:10.5f}{1e3 * (g - apse.dv):10.4f}")
        report(abs(g - apse.dv) < 1e-6, "apse-to-apse apogee burn within 1 mm/s")

    print("\n=== 2. GMAT flies the elliptical burn plans; arrival vs target orbit ===")
    for m in ELLIPTICAL:
        for o in evaluate(m):
            if not o.applicable:
                continue
            yf = gmat(console, lambda rep: fly_script(o.y0, o.vec_plan, rep))
            if yf is None:
                report(False, f"{o.name}: GMAT run failed")
                continue
            dh, de = target_mismatch(yf, m.r2, m.ra2, math.radians(m.di_deg), math.radians(m.dw_deg))
            report(dh < 1e-6 and de < 1e-6,
                   f"{m.alt1:.0f}x{m.ra1 - RE:.0f} -> {m.alt2:.0f}x{m.ra2 - RE:.0f} km di={m.di_deg:g} dw={m.dw_deg:g}: "
                   f"{o.name:32s} |dh|/h = {dh:.1e}, |de| = {de:.1e}")
    print(f"\n{'ALL CHECKS PASSED' if not fails else f'{fails} CHECK(S) FAILED'}")
    sys.exit(fails)


if __name__ == "__main__":
    main()
