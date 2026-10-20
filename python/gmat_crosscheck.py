"""Cross-check against NASA GMAT: GMAT's differential corrector solves the same
Hohmann transfers that odst.py sizes analytically, then the two are compared.

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

from odst import MU, RE, hohmann, vcirc, vvis

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
TEMPLATE = Template((HERE.parent / "gmat" / "hohmann_plane_change.script").read_text())
POINT_MASS = "FM.PointMasses = {Earth};"
J2 = ("FM.PrimaryBodies = {Earth};\nFM.GravityField.Earth.Degree = 2;\n"
      "FM.GravityField.Earth.Order = 0;\nFM.GravityField.Earth.PotentialFile = 'JGM2.cof';")
INC1 = 30.0   # start inclination; target is INC1 - di (see the script's header)


def run_gmat(console, alt1, alt2, di_deg, gravity, sign):
    r1, r2, di = RE + alt1, RE + alt2, math.radians(di_deg)
    vt1, vt2 = vvis(r1, (r1 + r2) / 2), vvis(r2, (r1 + r2) / 2)
    with tempfile.TemporaryDirectory() as tmp:
        report = Path(tmp) / "report.txt"
        script = Path(tmp) / "case.script"
        script.write_text(TEMPLATE.substitute(
            r1=f"{r1:.6f}", r2=f"{r2:.6f}", inc1=INC1, inc2=INC1 - di_deg, gravity=gravity,
            dv1=f"{vt1 - vcirc(r1):.6f}", dv2v=f"{vcirc(r2) * math.cos(di) - vt2:.6f}",
            dv2n=f"{sign * vcirc(r2) * math.sin(di):.6f}", report=report.as_posix()))
        subprocess.run([console, str(script)], cwd=Path(console).parent,
                       capture_output=True, text=True, timeout=600)
        if not report.exists():
            return None
        vals = [float(x) for x in report.read_text().split()]
    toi, goi, t, sma, ecc, inc = vals[0:3], vals[3:6], vals[6], vals[7], vals[8], vals[9]
    return {"dv": math.hypot(*toi) + math.hypot(*goi), "tof": t, "sma": sma, "ecc": ecc, "inc": inc,
            "inc_target": INC1 - di_deg}


def solve(console, alt1, alt2, di_deg, gravity):
    """The normal component's sign depends on the node; try both, keep the one that converged."""
    for sign in (1, -1):
        g = run_gmat(console, alt1, alt2, di_deg, gravity, sign)
        if g and abs(g["inc"] - g["inc_target"]) < 1e-6 and abs(g["sma"] - (RE + alt2)) < 1e-4:
            return g
    return None


def main():
    console = (sys.argv[1] if len(sys.argv) > 1 else
               os.environ.get("GMAT_CONSOLE", "D:/GMAT/bin/GmatConsole.exe"))
    cases = [("400 km -> GEO, di 28.5 deg, two-body", 400, 35786, 28.5, POINT_MASS),
             ("400 km -> GEO, coplanar, two-body", 400, 35786, 0.0, POINT_MASS),
             ("400 km -> GEO, di 28.5 deg, J2 (JGM-2)", 400, 35786, 28.5, J2)]
    fails = 0
    print(f"GMAT: {console}\n")
    print(f"{'Case':40s}{'ODST dv':>10s}{'GMAT dv':>10s}{'diff m/s':>10s}{'ODST TOF':>10s}{'GMAT TOF':>10s}{'diff s':>9s}")
    for name, alt1, alt2, di, grav in cases:
        dv, tof = hohmann(RE + alt1, RE + alt2, math.radians(di))[:2]
        g = solve(console, alt1, alt2, di, grav)
        if g is None:
            fails += 1
            print(f"{name:40s}  GMAT did not converge")
            continue
        ddv, dt = 1e3 * (g["dv"] - dv), g["tof"] - tof
        print(f"{name:40s}{dv:10.5f}{g['dv']:10.5f}{ddv:10.4f}{tof / 3600:9.4f}h{g['tof'] / 3600:9.4f}h{dt:9.3f}")
        if grav == POINT_MASS:        # same physics: must agree to 1 mm/s and 0.01 s
            ok = abs(ddv) < 1e-3 and abs(dt) < 0.01
            fails += not ok
            print(f"  [{'PASS' if ok else 'FAIL'}] two-body agreement within 1 mm/s and 0.01 s"
                  f" (GMAT final e = {g['ecc']:.1e})")
    print(f"\n{'ALL CHECKS PASSED' if not fails else f'{fails} CHECK(S) FAILED'}")
    sys.exit(fails)


if __name__ == "__main__":
    main()
