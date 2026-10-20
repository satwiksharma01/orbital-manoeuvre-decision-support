"""Validation of the ODST models: textbook benchmarks + numerical propagation."""
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from odst import (G0, MU, RE, Mission, apse_transfer, bielliptic, evaluate, gravity_loss, hohmann,
                  hohmann_optimal_split, lambert, lambert_optimal, one_tangent, orbit_elements,
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
       Mission(alt1=500, apo1=5000, alt2=800, apo2=12000, dw_deg=60, di_deg=10)]
r = lambda a: RE + a
coarse = lambert_optimal(r(500), r(5000), r(800), r(12000), math.radians(10), math.radians(60), 1e9)["dv"]
fine = lambert_optimal(r(500), r(5000), r(800), r(12000), math.radians(10), math.radians(60), 1e9, n_grid=48)["dv"]
check("Lambert search, 16x16 grid vs 48x48 grid, misaligned + 10 deg [m/s]", 1e3 * coarse, 1e3 * fine, 0.01)
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

print("\n=== 6. Mission case studies ===")
for m in cases + ell:
    print_table(m, evaluate(m))

print(f"\n{'ALL CHECKS PASSED' if fails == 0 else f'{fails} CHECK(S) FAILED'}")
sys.exit(fails)
