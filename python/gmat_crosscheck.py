"""Cross-check against NASA GMAT.

1. GMAT's differential corrector solves the same Hohmann transfers that odst.py sizes
   analytically, and the GTO -> GEO apogee burn of the apse-to-apse model.
2. GMAT runs the scripts the tool exports (odst.gmat_script) for every impulsive option of
   circular and elliptical missions, some given by orbital elements, and the orbit it arrives
   on is compared with the target: its angular-momentum and eccentricity vectors in
   EarthMJ2000Eq, so the plane and the apse line count, not only the shape.
3. With the spacecraft model the script carries the spacecraft's mass and a chemical tank, and
   GMAT draws each burn's propellant from it: GMAT's propellant used vs the tool's.
4. Uncertainty: GMAT flies each burn off course by the pointing error of a Monte Carlo run,
   then the trim, and must arrive on target having used the propellant the tool's trim cost,
   2 dv sin(theta/2) per burn, predicts.
5. J2: GMAT flies the tool's J2 plans (the nodal drift, and transfers aimed at where the target's
   plane will be) in a J2 field with the target orbit alongside, and compares the two planes.

Needs a local GMAT (https://sourceforge.net/projects/gmat/). Not run in CI.
Run: python python/gmat_crosscheck.py [path/to/GmatConsole.exe]
     (default: $GMAT_CONSOLE, else D:/GMAT/bin/GmatConsole.exe)
"""
import math
import os
import subprocess
import sys
import tempfile
from dataclasses import replace
from pathlib import Path
from string import Template

import numpy as np

from odst import (G0, MU, RE, Mission, _rotate, evaluate, gmat_script, hohmann, mc_normals, orbit_basis, perifocal, vcirc,
                  vvis)

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
TEMPLATE = Template((HERE.parent / "gmat" / "hohmann_plane_change.script").read_text())
APOGEE = Template((HERE.parent / "gmat" / "apogee_burn.script").read_text())
POINT_MASS = "FM.PointMasses = {Earth};"
J2 = ("FM.PrimaryBodies = {Earth};\nFM.GravityField.Earth.Degree = 2;\n"
      "FM.GravityField.Earth.Order = 0;\nFM.GravityField.Earth.PotentialFile = 'JGM2.cof';")
INC1 = 30.0   # start inclination; target is INC1 - di (see the scripts' headers)
CIRCULAR = [Mission(),                                                                   # case A: LEO -> GEO
            Mission(alt1=7000 - RE, alt2=105000 - RE, di_deg=0, tof_max_h=24 * 60),       # case B: R = 15
            Mission(alt1=400, alt2=800, di_deg=0, tof_max_h=0.6),                        # case C: 36 min
            Mission(alt1=400, alt2=400, di_deg=60),                                      # plane change only
            Mission(alt1=35786, alt2=400, di_deg=28.5, tof_max_h=3)]                     # lowering, fast
ELLIPTICAL = [Mission(alt1=250, apo1=35786, alt2=35786, di_deg=28.5),                   # GTO -> GEO
              Mission(alt1=400, alt2=600, apo2=39700, dw_deg=270, di_deg=0),            # LEO -> Molniya
              Mission(alt1=500, apo1=5000, alt2=800, apo2=12000, dw_deg=60, di_deg=0),  # apse lines 60 deg apart
              Mission(alt1=500, apo1=5000, alt2=800, apo2=12000, dw_deg=60, di_deg=10)]
ELEMENTS = [Mission(alt1=400, alt2=35786, elements=(51.6, 120, 0, 0, 0, 0)),                   # 51.6 deg LEO -> GEO
            Mission(alt1=500, alt2=500, elements=(51.6, 0, 0, 51.6, 10, 0)),                   # RAAN change only
            Mission(alt1=250, apo1=35786, alt2=35786, elements=(6, 0, 178, 0, 0, 0)),          # GTO, apogee off the node
            Mission(alt1=400, alt2=600, apo2=39700, elements=(51.6, 30, 0, 63.4, 100, 270), tof_max_h=30)]

SPACECRAFT = [Mission(mass0_kg=1000),                                                       # LEO -> GEO, defaults
              Mission(alt1=250, apo1=35786, alt2=35786, di_deg=28.5, mass0_kg=2000, thrust_chem_N=2000,
                      prop_chem_kg=1200, isp_chem=310)]                                     # GTO -> GEO


def target_vectors(m):
    """The target's unit vectors to periapsis and along its normal in EarthMJ2000Eq."""
    if m.elements is not None:
        P, _, W = perifocal(*map(math.radians, m.elements[3:]))
        return P, W
    P, Q = orbit_basis(math.radians(m.di_deg), math.radians(m.dw_deg))
    return P, np.cross(P, Q)


def gmat(console, make_script):
    """Run the script make_script(report_path) builds; return the reported numbers or None."""
    with tempfile.TemporaryDirectory() as tmp:
        report, script = Path(tmp) / "report.txt", Path(tmp) / "case.script"
        script.write_text(make_script(report.as_posix()))
        subprocess.run([console, str(script)], cwd=Path(console).parent,
                       capture_output=True, text=True, timeout=600)
        return [float(x) for x in report.read_text().splitlines()[-1].split()] if report.exists() else None


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

    print("\n=== 2. GMAT runs the exported scripts; arrival vs target orbit ===")
    for m in CIRCULAR + ELLIPTICAL + ELEMENTS:
        P, W = target_vectors(m)
        p, ecc = 2 * m.r2 * m.ra2 / (m.r2 + m.ra2), (m.ra2 - m.r2) / (m.ra2 + m.r2)
        for o in evaluate(m):
            if not o.vec_plan or not o.applicable:
                continue
            g = gmat(console, lambda rep: gmat_script(m, o, rep))
            orbits = (f"elements {m.elements}" if m.elements else
                      f"{m.alt1:.0f}x{m.ra1 - RE:.0f} -> {m.alt2:.0f}x{m.ra2 - RE:.0f} km di={m.di_deg:g} dw={m.dw_deg:g}")
            name = f"{orbits}: {o.name:36s}"
            if g is None:
                report(False, f"{name} GMAT run failed")
                continue
            r, v = np.array(g[:3]), np.array(g[3:6])
            h = np.cross(r, v)
            dh = np.linalg.norm(h - math.sqrt(MU * p) * W) / math.sqrt(MU * p)
            de = np.linalg.norm(np.cross(v, h) / MU - r / np.linalg.norm(r) - ecc * P)
            gm = f"; GMAT i {g[8]:.6f}, RAAN {g[9]:.6f}" if m.elements else ""
            report(dh < 1e-6 and de < 1e-6, f"{name} |dh|/h = {dh:.1e}, |de| = {de:.1e}{gm}")

    print("\n=== 3. GMAT draws the propellant from a tank (spacecraft model) ===")
    for m in SPACECRAFT:
        for o in evaluate(m):
            if not o.vec_plan or not o.applicable:
                continue
            g = gmat(console, lambda rep: gmat_script(m, o, rep))
            if g is None:
                report(False, f"{o.name}: GMAT run failed")
                continue
            used = m.prop_chem_kg - g[11]
            ref = m.mass0_kg * (1 - math.exp(-o.dv / (m.isp_chem * G0)))     # impulsive: no gravity loss, no margin
            report(abs(used - ref) < 1e-6, f"{m.mass0_kg:g} kg, {m.alt1:.0f}x{m.ra1 - RE:.0f} -> {m.alt2:.0f} km: {o.name:36s} "
                                           f"propellant used, GMAT {used:.6f} kg vs tool {ref:.6f} kg")

    print("\n=== 4. Uncertainty: off-course burns of three Monte Carlo runs (pointing 2 deg, 1 sigma), each trimmed ===")
    Z, sig = mc_normals(), math.radians(2.0)
    for m in SPACECRAFT:
        P, W = target_vectors(m)
        p, ecc = 2 * m.r2 * m.ra2 / (m.r2 + m.ra2), (m.ra2 - m.r2) / (m.ra2 + m.r2)
        for o in evaluate(m):
            if not o.vec_plan or not o.applicable:
                continue
            for run in range(3):
                plan, trims = [], 0.0
                for j, (t, dv) in enumerate(o.vec_plan):        # turn the burn by theta about an axis across it
                    z1, z2 = Z[run, 2 + 2 * j], Z[run, 3 + 2 * j]
                    e1 = np.cross(dv, [0.0, 0.0, 1.0] if abs(dv[2]) < 0.9 * np.linalg.norm(dv) else [1.0, 0.0, 0.0])
                    e1 /= np.linalg.norm(e1)
                    e2 = np.cross(dv / np.linalg.norm(dv), e1)
                    flown = _rotate(dv, z1 * e1 + z2 * e2, sig * math.hypot(z1, z2))
                    plan += [(t, flown), (t, dv - flown)]
                    trims += 2 * np.linalg.norm(dv) * math.sin(sig * math.hypot(z1, z2) / 2)
                g = gmat(console, lambda rep: gmat_script(m, replace(o, vec_plan=plan), rep))
                if g is None:
                    report(False, f"{o.name}, run {run + 1}: GMAT run failed")
                    continue
                r, v = np.array(g[:3]), np.array(g[3:6])
                h = np.cross(r, v)
                dh = np.linalg.norm(h - math.sqrt(MU * p) * W) / math.sqrt(MU * p)
                de = np.linalg.norm(np.cross(v, h) / MU - r / np.linalg.norm(r) - ecc * P)
                used = m.prop_chem_kg - g[11]
                ref = m.mass0_kg * (1 - math.exp(-(o.dv + trims) / (m.isp_chem * G0)))
                report(dh < 1e-6 and de < 1e-6 and abs(used - ref) < 1e-6,
                       f"{m.mass0_kg:g} kg, {m.alt1:.0f}x{m.ra1 - RE:.0f} -> {m.alt2:.0f} km: {o.name:36s} run {run + 1}: trims "
                       f"{1e3 * trims:6.1f} m/s; |dh|/h {dh:.1e}, |de| {de:.1e}; propellant GMAT {used:.6f} vs tool {ref:.6f} kg")
    print("\n=== 5. J2: the tool's J2 plans in GMAT's J2 field (JGM-2), with the target orbit alongside ===")

    def planes(g, k):
        """Angle (deg) between Sat's plane (INC, RAAN at g[8], g[9]) and TargetSat's (at g[k], g[k + 1])."""
        n = lambda i, O: np.array([math.sin(O) * math.sin(i), -math.cos(O) * math.sin(i), math.cos(i)])
        return math.degrees(math.acos(min(1.0, n(*map(math.radians, g[8:10])) @ n(*map(math.radians, g[k:k + 2])))))
    for kw in (dict(alt1=500, alt2=500, elements=(51.6, 0, 0, 51.6, 10, 0)), dict(alt1=500, alt2=500, elements=(51.6, 0, 0, 53, 10, 0)),
               dict(alt1=400, alt2=1200, elements=(51.6, 15, 0, 51.6, 0, 0))):
        m = Mission(j2=True, tof_max_h=16 * 24, **kw)
        o = next(o for o in evaluate(m) if o.key == "nodal_drift")
        g = gmat(console, lambda rep: gmat_script(m, o, rep))
        if g is None:
            report(False, f"nodal drift {kw['elements']}: GMAT run failed")
            continue
        ang = planes(g, 11)
        report(ang < 0.02 * m.di_deg, f"Nodal drift {m.alt1:.0f} -> {m.alt2:.0f} km, i {m.elements[0]:g} -> {m.elements[3]:g}, RAAN "
                                      f"{m.elements[1]:g} -> {m.elements[4]:g}: {o.dv:.4f} km/s, wait {o.params['wait_days']:.2f} days at "
                                      f"{o.params['alt_d']:.0f} km; GMAT: {ang:.3f} deg from the target's plane (direct plane change "
                                      f"{m.di_deg:.2f} deg), SMA {g[6]:.1f} vs {g[13]:.1f} km")
    for kw in (dict(alt1=20200, alt2=500, elements=(55, 0, 0, 55, 3, 0)), dict(alt1=35786, alt2=600, elements=(30, 0, 0, 51.6, 20, 0)),
               dict(alt1=400, alt2=20200, elements=(28.5, 0, 0, 55, 40, 0))):
        m = Mission(j2=True, tof_max_h=1e5, **kw)
        o = next(o for o in evaluate(m) if o.key == "hohmann")
        o2 = next(o for o in evaluate(replace(m, j2=False)) if o.key == "hohmann")
        ga, gp = (gmat(console, lambda rep, oo=oo: gmat_script(m, oo, rep)) for oo in (o, o2))
        if ga is None or gp is None:
            report(False, f"Hohmann {kw['elements']}: GMAT run failed")
            continue
        aimed, plain = planes(ga, 11), planes(gp, 11)
        descent = m.r2 < m.r1      # the target's own drift dominates: the aim must win clearly; 0.1 deg: J2 within one orbit, of order J2 (R/p)^2
        report(aimed < 0.1 and aimed < plain / 4 if descent else aimed < 0.15,
               f"Hohmann {m.alt1:.0f} -> {m.alt2:.0f} km, i {m.elements[0]:g} -> {m.elements[3]:g}, {o.tof / 3600:.1f} h: GMAT, two-body plan "
               f"{plain:.3f} deg from the target's plane, aimed {aimed:.4f} deg ({1e3 * o.j2_dv:+.1f} m/s)")
    print(f"\n{'ALL CHECKS PASSED' if not fails else f'{fails} CHECK(S) FAILED'}")
    sys.exit(fails)


if __name__ == "__main__":
    main()
