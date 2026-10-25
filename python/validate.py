"""Validation of the ODST models: textbook benchmarks + numerical propagation."""
import json
import math
import re
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from odst import (G0, H_ARC_MIN, MU, RE, Mission, apse_transfer, arc_min_radius, bielliptic, burns_json, evaluate,
                  gmat_script, gravity_loss, hohmann, hohmann_optimal_split, inertial_plan, kepler, lambert,
                  lambert_optimal, low_thrust, one_tangent, orbit_basis, orbit_elements, orbit_state, perifocal,
                  print_table, simulate, simulate_vec, target_mismatch, vcirc, vvis)

sys.stdout.reconfigure(encoding="utf-8")
fails = 0


def check(name, got, ref, tol):
    global fails
    ok = abs(got - ref) <= tol
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: got {got:.5f}, ref {ref:.5f} (tol {tol})")


print("=== 1. Textbook benchmarks ===")
# Standard textbook benchmark: 300 km LEO -> GEO coplanar Hohmann
dv, tof, _, _ = hohmann(RE + 300, 42164.0, 0.0)
check("Hohmann 300 km->GEO dv [km/s]", dv, 3.893, 0.005)
check("Hohmann 300 km->GEO TOF [h]", tof / 3600, 5.27, 0.02)
# Bi-elliptic/Hohmann crossover: R = 11.94 equal dv as rB -> infinity
r1 = 7000.0
for R, expect_bi_better in [(11.0, False), (15.0, True), (16.0, True)]:
    r2 = R * r1
    dh = hohmann(r1, r2, 0)[0]
    db = bielliptic(r1, r2, 1e4 * r2, 0)[0]
    ok = (db < dh) == expect_bi_better
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] R={R}: Hohmann {dh:.4f} vs bi-elliptic(rB->inf) {db:.4f}")
# Pure plane change: dv = 2 v sin(di/2)
dv = hohmann(7000, 7000, math.radians(30))[0]
check("Simple plane change 30 deg at 7000 km", dv, 2 * vcirc(7000) * math.sin(math.radians(15)), 1e-9)
# One-tangent with k -> 1 must reduce to Hohmann
check("One-tangent(k=1) == Hohmann dv", one_tangent(7000, 20000, 1 + 1e-10, 0)[0],
      hohmann(7000, 20000, 0)[0], 1e-4)

print("\n=== 2. Numerical verification of burn plans (two-body propagation) ===")
cases = [
    Mission(alt1=400, alt2=35786, di_deg=28.5),
    Mission(alt1=621.863, alt2=15 * (RE + 621.863) - RE, di_deg=0, tof_max_h=24 * 60),
    Mission(alt1=400, alt2=800, di_deg=0, tof_max_h=0.6),
    Mission(alt1=35786, alt2=400, di_deg=10),   # lowering
]
for m in cases:
    for o in evaluate(m):
        if not o.burn_plan or not o.applicable:
            continue
        _, _, yf, dvsim = simulate(m.r1, o.burn_plan, coast_after=3600)
        a, e, inc = orbit_elements(yf)
        good = abs(a - m.r2) < 0.5 and e < 1e-4 and abs(inc - m.di_deg) < 1e-3 and abs(dvsim - o.dv) < 1e-4
        fails += not good
        print(f"[{'PASS' if good else 'FAIL'}] {m.alt1:.0f}->{m.alt2:.0f} km di={m.di_deg}: {o.name:30s} "
              f"a={a:10.2f} (target {m.r2:.2f}) e={e:.1e} i={inc:7.3f} dv_sim={dvsim:.5f} dv_model={o.dv:.5f}")

print("\n=== 3. Δv margin (low thrust needs 5.898 km/s for LEO -> GEO) ===")
for budget, margin, expect in [(6.0, 0, True), (6.0, 5, False), (6.5, 5, True)]:
    lt = next(o for o in evaluate(Mission(dv_budget=budget, margin_pct=margin)) if o.key == "low_thrust")
    ok = lt.feasible == expect
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] budget {budget} km/s, margin {margin}%: low thrust feasible = {lt.feasible}")

print("\n=== 4. Gravity loss: Robbins estimate vs numerical finite burn ===")


def finite_apoapsis(r1, tau, tw, isp):
    """Constant thrust at a fixed attitude (along the impulsive burn), centred on the burn point."""
    ve, a0 = isp * G0, tw * G0
    n, v = math.sqrt(MU / r1 ** 3), vcirc(r1)
    th = -n * tau / 2
    y0 = [r1 * math.cos(th), r1 * math.sin(th), -v * math.sin(th), v * math.cos(th)]
    f = lambda t, y: [y[2], y[3], -MU * y[0] / math.hypot(y[0], y[1]) ** 3,
                      -MU * y[1] / math.hypot(y[0], y[1]) ** 3 + a0 / (1 - a0 / ve * t)]
    y = solve_ivp(f, (0, tau), y0, rtol=1e-11, atol=1e-11).y[:, -1]
    r, eps = math.hypot(y[0], y[1]), (y[2] ** 2 + y[3] ** 2) / 2 - MU / math.hypot(y[0], y[1])
    if eps >= 0:
        return math.inf                                    # escaped: overshoot
    a = -MU / (2 * eps)
    h = y[0] * y[3] - y[1] * y[2]
    return a * (1 + math.sqrt(max(0.0, 1 - h * h / (MU * a))))


r1, r2, isp = RE + 400, 42164.0, 320.0
dv1 = vvis(r1, (r1 + r2) / 2) - vcirc(r1)
for tw in (0.5, 0.2):
    est = gravity_loss([(r1, dv1)], isp, tw)
    tau0 = isp / tw * (1 - math.exp(-dv1 / (isp * G0)))
    tau = brentq(lambda t: finite_apoapsis(r1, t, tw, isp) - r2, tau0, 1.3 * tau0, xtol=1e-6)
    num = -isp * G0 * math.log(1 - tw / isp * tau) - dv1
    check(f"LEO->GEO first burn loss, T/W {tw} [m/s] (Robbins vs numerical)", 1e3 * est, 1e3 * num, 0.1 * 1e3 * num)
o = next(o for o in evaluate(Mission(tw_chem=0.2)) if o.key == "hohmann")
check("Propellant fraction includes the gravity loss (Hohmann LEO->GEO, T/W 0.2)", o.prop_frac,
      1 - math.exp(-(o.dv + o.dv_loss) / (o.isp * G0)), 1e-12)

print("\n=== 4b. Edelbaum low thrust vs the averaged equations, integrated with Edelbaum's yaw law ===")
for r1, r2, di_deg in ((RE + 400, 42164.0, 28.5), (42164.0, RE + 400, 10.0), (7000.0, 7000.0, 30.0)):
    f, di = 1e-6, math.radians(di_deg)                    # 1e-6 km/s^2 = 1 mm/s^2
    dv, tof = low_thrust(r1, r2, di, f * 1e3)
    V0, V1 = vcirc(r1), vcirc(r2)
    b0 = math.atan2(math.sin(math.pi / 2 * di), V0 / V1 - math.cos(math.pi / 2 * di))
    beta = lambda t: math.atan2(V0 * math.sin(b0), V0 * math.cos(b0) - f * t)
    end = solve_ivp(lambda t, y: [-f * math.cos(beta(t)), 2 / math.pi * f * math.sin(beta(t)) / y[0]], (0, tof),
                    [V0, 0.0], rtol=1e-11, atol=1e-13).y[:, -1]
    check(f"Edelbaum {r1:.0f} -> {r2:.0f} km, di {di_deg}: speed after the closed-form TOF [km/s]", end[0], V1, 1e-8)
    check(f"Edelbaum {r1:.0f} -> {r2:.0f} km, di {di_deg}: plane change after the closed-form TOF [deg]", math.degrees(end[1]), di_deg, 1e-6)
lt = next(o for o in evaluate(Mission(di_deg=120)) if o.key == "low_thrust")
ok = not lt.applicable
fails += not ok
print(f"[{'PASS' if ok else 'FAIL'}] Edelbaum not offered beyond 114.6 deg (di 120): applicable = {lt.applicable} ({lt.note})")

print("\n=== 5. Elliptical orbits and Lambert ===")
# Curtis, Orbital Mechanics for Engineering Students, Example 5.2
r1v, r2v = np.array([5000.0, 10000.0, 2100.0]), np.array([-14600.0, 2500.0, 7000.0])
v1, v2 = lambert(r1v, r2v, 3600.0)
for k, (got, ref) in enumerate(zip(np.concatenate([v1, v2]), [-5.9925, 1.9254, 3.2456, -3.3125, -4.1966, -0.38529])):
    check(f"Lambert, Curtis Ex. 5.2, {'v1' if k < 3 else 'v2'}[{k % 3}] [km/s]", got, ref, 1e-4)
yf = solve_ivp(lambda t, y: np.concatenate([y[3:], -MU * y[:3] / np.linalg.norm(y[:3]) ** 3]), (0, 3600.0),
               np.concatenate([r1v, v1]), rtol=1e-12, atol=1e-12).y[:, -1]
check("Lambert v1 propagated for 1 h reaches r2 [km]", np.linalg.norm(yf[:3] - r2v), 0.0, 1e-6)
di = math.radians(20)
check("Apse-to-apse between circular orbits == Hohmann, optimal split [km/s]",
      apse_transfer(7000, 7000, 14000, 14000, di, 0)[0]["dv"],
      hohmann(7000, 14000, di, hohmann_optimal_split(7000, 14000, di))[0], 1e-9)
check("Lambert optimum, coplanar 7000 -> 14000 km, vs Hohmann [m/s]",
      1e3 * lambert_optimal(7000, 7000, 14000, 14000, 0, 0, 1e9)["dv"], 1e3 * hohmann(7000, 14000, 0)[0], 1.0)
ell = [Mission(alt1=250, apo1=35786, alt2=35786, di_deg=28.5),                  # GTO -> GEO
       Mission(alt1=400, alt2=600, apo2=39700, dw_deg=270, di_deg=0),           # LEO -> Molniya
       Mission(alt1=500, apo1=5000, alt2=800, apo2=12000, dw_deg=60, di_deg=0), # apse lines 60 deg apart
       Mission(alt1=500, apo1=5000, alt2=800, apo2=12000, dw_deg=60, di_deg=10),
       Mission(alt1=300, apo1=12500, alt2=350, apo2=8600, dw_deg=90, di_deg=0, tof_max_h=1950 / 3600)]  # 100 km floor active
r = lambda a: RE + a
coarse = lambert_optimal(r(500), r(5000), r(800), r(12000), math.radians(10), math.radians(60), 1e9)["dv"]
fine = lambert_optimal(r(500), r(5000), r(800), r(12000), math.radians(10), math.radians(60), 1e9, n_grid=48)["dv"]
check("Lambert search, 16x16 grid vs 48x48 grid, misaligned + 10 deg [m/s]", 1e3 * coarse, 1e3 * fine, 0.01)
# Lowest point of a conic arc: closed form vs a propagated trajectory, refined to the true minimum
rng, worst = np.random.default_rng(3), 0.0
for _ in range(40):
    a, b = rng.uniform(6700, 30000, 2)
    u1, u2 = (v / np.linalg.norm(v) for v in rng.normal(size=(2, 3)))
    dt = math.sqrt((a + b) ** 3 / MU) * rng.uniform(0.05, 2.5)
    s = lambert(a * u1, b * u2, dt)
    if s is None:
        continue
    run = solve_ivp(lambda t, y: np.concatenate([y[3:], -MU * y[:3] / np.linalg.norm(y[:3]) ** 3]), (0, dt),
                    np.concatenate([a * u1, s[0]]), rtol=1e-12, atol=1e-12, dense_output=True)
    ts = np.linspace(0, dt, 2001)
    k = int(np.argmin([np.linalg.norm(run.sol(t)[:3]) for t in ts]))
    num = minimize_scalar(lambda t: np.linalg.norm(run.sol(t)[:3]), bounds=(ts[max(k - 1, 0)], ts[min(k + 1, 2000)]),
                          method="bounded", options={"xatol": 1e-9}).fun
    num = min(num, a, b)                     # the bounded search never lands exactly on an end point
    worst = max(worst, abs(num - arc_min_radius(a * u1, s[0], b * u2)))
check("Lowest point of 40 random transfer arcs, closed form vs propagation [km]", worst, 0.0, 1e-6)
# 100 km floor: without it this mission's optimum dives 327 km below the surface (audit, 2026-10-22).
# Global optimum with the floor, from differential evolution (8 runs agree): 2.546787011 km/s.
m = ell[-1]
o = next(o for o in evaluate(m) if o.key == "lambert")
y = np.array(o.y0, float); y[3:] += o.vec_plan[0][1]
_, _, y2, _ = simulate_vec(y, [], coast_after=o.vec_plan[1][0])
check("Floor-active Lambert: the arc's lowest point sits on the 100 km floor [km altitude]",
      arc_min_radius(y[:3], y[3:], y2[:3]) - RE, H_ARC_MIN, 1e-3)
check("Floor-active Lambert vs the global optimum [m/s]", 1e3 * o.dv, 2546.787011, 0.01)
for m in ell:
    for o in evaluate(m):
        if not o.applicable:
            continue
        _, _, yf, dvsim = simulate_vec(o.y0, o.vec_plan, coast_after=3600)
        dh, de = target_mismatch(yf, m.r2, m.ra2, math.radians(m.di_deg), math.radians(m.dw_deg))
        good = dh < 1e-6 and de < 1e-6 and abs(dvsim - o.dv) < 1e-9
        fails += not good
        print(f"[{'PASS' if good else 'FAIL'}] {m.alt1:.0f}x{m.ra1 - RE:.0f} -> {m.alt2:.0f}x{m.ra2 - RE:.0f} km "
              f"di={m.di_deg} dw={m.dw_deg}: {o.name:32s} dv={o.dv:.6f} TOF={o.tof / 3600:7.3f} h  "
              f"arrival |dh|/h={dh:.1e} |de|={de:.1e}")

print("\n=== 5b. Export: exact burn plans, GMAT script and JSON ===")
# Kepler propagation vs Kepler's equation: start and end points from eccentric anomalies
rng = np.random.default_rng(7)
worst = 0.0
for _ in range(200):
    rp = RE + rng.uniform(200, 2000)
    ra = rp * rng.uniform(1, 60)
    a, e = (rp + ra) / 2, (ra - rp) / (ra + rp)
    P, Q = orbit_basis(rng.uniform(0, math.pi), rng.uniform(0, 2 * math.pi))
    th = lambda E: 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))
    t = lambda E: math.sqrt(a ** 3 / MU) * (E - e * math.sin(E))
    E0, E1 = rng.uniform(0, 2 * math.pi, 2)
    E1 += 2 * math.pi * (E1 < E0)
    y = kepler(np.concatenate(orbit_state(rp, ra, P, Q, th(E0))), t(E1) - t(E0))
    r, v = orbit_state(rp, ra, P, Q, th(E1))
    worst = max(worst, np.linalg.norm(y[:3] - r) / np.linalg.norm(r), np.linalg.norm(y[3:] - v) / np.linalg.norm(v))
check("Kepler propagation vs Kepler's equation, 200 random arcs up to e = 0.97 [1e-12 relative]", 1e12 * worst, 0.0, 100)
# The circular models' exact vector plans: propagated, they reach the target plane itself (normal
# (0, -sin di, cos di)), not only its inclination, at the closed-form cost
for m in cases + [Mission(alt1=400, alt2=400, di_deg=60), Mission(alt1=35786, alt2=400, di_deg=28.5, tof_max_h=3)]:
    for o in evaluate(m):
        if not o.burn_plan or not o.applicable:
            continue
        _, _, yf, dvsim = simulate_vec(o.y0, o.vec_plan, coast_after=3600)
        dh, de = target_mismatch(yf, m.r2, m.ra2, math.radians(m.di_deg), 0.0)
        good = dh < 1e-6 and de < 1e-6 and abs(dvsim - o.dv) < 1e-12
        fails += not good
        print(f"[{'PASS' if good else 'FAIL'}] {m.alt1:.0f}->{m.alt2:.0f} km di={m.di_deg}: {o.name:30s} vector plan "
              f"|dh|/h={dh:.1e} |de|={de:.1e}, sum|dv| - dv={dvsim - o.dv:.1e}")
# GMAT script and JSON read back to the exact plan (GMAT flies them: gmat_crosscheck.py)
for m in (cases[0], ell[1], Mission(alt1=400, alt2=600, apo2=39700, elements=(51.6, 30, 0, 63.4, 100, 270), tof_max_h=30)):
    for o in evaluate(m):
        if not o.vec_plan:
            continue
        txt, js = gmat_script(m, o), json.loads(burns_json(m, o))
        y0i, plani = inertial_plan(m, o)
        val = lambda k: float(re.search(rf"^{re.escape(k)} = (\S+);$", txt, re.M).group(1))
        plan, t = [], 0.0
        for line in txt[txt.index("BeginMissionSequence"):].splitlines():
            if line.startswith("Propagate") and len(plan) < len(plani):
                t += float(re.search(r"ElapsedSecs = (\S+)}", line).group(1))
            elif line.startswith("Maneuver"):
                k = len(plan) + 1
                plan.append((t, [val(f"Burn{k}.Element{j}") for j in (1, 2, 3)]))
        y0 = [val(f"Sat.{k}") for k in ("X", "Y", "Z", "VX", "VY", "VZ")]
        jplan = [(b["t_s"], b["dv_kms"]) for b in js["burns"]]
        good = (val("Earth.Mu") == MU and y0 == list(y0i) and js["start_state_km_kms"] == list(y0i)
                and all(abs(p[0] - q[0]) < 1e-9 * max(1.0, q[0]) and p[1] == list(q[1]) for p, q in zip(plan, plani))
                and jplan == [(q[0], list(q[1])) for q in plani] and len(plan) == len(plani)
                and js["elements_deg"] == (None if m.elements is None else list(m.elements)))
        fails += not good
        print(f"[{'PASS' if good else 'FAIL'}] {o.name:32s} GMAT script and JSON give back the plan "
              f"({len(plan)} burns, start state, Earth mu{', elements' if m.elements else ''})")

print("\n=== 5c. Orbital elements (inclination, RAAN, argument of periapsis) ===")
dvs = lambda m: {o.key: o.dv for o in evaluate(m) if o.applicable}
diff = lambda a, b: max(abs(a[k] - b[k]) for k in a) if a.keys() == b.keys() else math.inf
# Elements that put the start orbit in x-y with the node on +x are the relative inputs themselves
for rel, el in [(Mission(alt1=400, alt2=600, apo2=39700, di_deg=20, dw_deg=270), (0, 0, 0, 20, 0, 270)),
                (Mission(alt1=250, apo1=35786, alt2=35786, di_deg=28.5), (0, 0, 0, 28.5, 0, 0))]:
    m = Mission(alt1=rel.alt1, apo1=rel.apo1, alt2=rel.alt2, apo2=rel.apo2, elements=el)
    check(f"Elements {el} == relative inputs, every option [km/s]", diff(dvs(m), dvs(rel)), 0.0, 1e-12)
# Turning both orbits together (both RAANs + 40 deg) changes nothing
c = dvs(Mission(alt1=300, apo1=12000, alt2=800, apo2=20000, elements=(28.5, 10, 30, 40, 70, 120), tof_max_h=20))
d = dvs(Mission(alt1=300, apo1=12000, alt2=800, apo2=20000, elements=(28.5, 50, 30, 40, 110, 120), tof_max_h=20))
check("Both orbits turned 40 deg about the pole: Lambert, 3D [km/s]", diff(c, d), 0.0, 1e-9)
# A RAAN change alone: the planes are 2 asin(sin i sin(dRAAN/2)) apart (spherical trigonometry)
m = Mission(alt1=500, alt2=500, elements=(51.6, 0, 0, 51.6, 10, 0))
ref = 2 * math.asin(math.sin(math.radians(51.6)) * math.sin(math.radians(5)))
check("RAAN 0 -> 10 deg at i = 51.6 deg: angle between the planes [deg]", m.di_deg, math.degrees(ref), 1e-12)
check("  its single-burn plane change, 2 v sin(di/2) [km/s]", dvs(m)["hohmann"], 2 * vcirc(m.r1) * math.sin(ref / 2), 1e-12)
# GTO with its apogee on the node line equals the relative GTO; 2 deg off it, only Lambert is exact
g = dvs(Mission(alt1=250, apo1=35786, alt2=35786, elements=(6, 0, 180, 0, 0, 0)))
check("GTO i 6, argp 180 -> GEO: apse-to-apse == relative di 6 [km/s]", g["apse"],
      dvs(Mission(alt1=250, apo1=35786, alt2=35786, di_deg=6))["apse"], 1e-12)
o = evaluate(Mission(alt1=250, apo1=35786, alt2=35786, elements=(6, 0, 178, 0, 0, 0)))
ok = not o[0].applicable and o[1].applicable and o[1].dv > g["apse"]
fails += not ok
print(f"[{'PASS' if ok else 'FAIL'}] GTO i 6, argp 178 -> GEO: apse-to-apse not offered ({o[0].note}); "
      f"Lambert {o[1].dv:.6f} km/s, {1e3 * (o[1].dv - g['apse']):.1f} m/s above the aligned case")
# Every plan, turned into EarthMJ2000Eq and propagated, arrives on the target given by its own elements
for m in (Mission(alt1=400, alt2=35786, elements=(51.6, 120, 0, 0, 0, 0)),
          Mission(alt1=400, alt2=35786, elements=(28.5, 0, 0, 0, 0, 0), tof_max_h=3),
          Mission(alt1=500, alt2=500, elements=(51.6, 0, 0, 51.6, 10, 0)),
          Mission(alt1=250, apo1=35786, alt2=35786, elements=(6, 0, 178, 0, 0, 0)),
          Mission(alt1=400, alt2=600, apo2=39700, elements=(51.6, 30, 0, 63.4, 100, 270), tof_max_h=30),
          Mission(alt1=300, apo1=12000, alt2=800, apo2=20000, elements=(28.5, 10, 30, 40, 70, 120), tof_max_h=20)):
    P2, _, W2 = perifocal(*map(math.radians, m.elements[3:]))
    for o in evaluate(m):
        if not o.vec_plan or not o.applicable:
            continue
        _, _, yf, dvsim = simulate_vec(*inertial_plan(m, o), coast_after=3600)
        r, v = yf[:3], yf[3:]
        h, p = np.cross(r, v), 2 * m.r2 * m.ra2 / (m.r2 + m.ra2)
        dh = np.linalg.norm(h - math.sqrt(MU * p) * W2) / math.sqrt(MU * p)
        de = np.linalg.norm(np.cross(v, h) / MU - r / np.linalg.norm(r) - (m.ra2 - m.r2) / (m.ra2 + m.r2) * P2)
        good = dh < 1e-6 and de < 1e-6 and abs(dvsim - o.dv) < 1e-9
        fails += not good
        print(f"[{'PASS' if good else 'FAIL'}] {str(m.elements):32s} {o.name:32s} dv={o.dv:.6f}  arrival vs the "
              f"target's elements |dh|/h={dh:.1e} |de|={de:.1e}")

print("\n=== 6. Mission case studies ===")
for m in cases + ell:
    print_table(m, evaluate(m))

print(f"\n{'ALL CHECKS PASSED' if fails == 0 else f'{fails} CHECK(S) FAILED'}")
sys.exit(fails)
