"""Validation of the ODST models: textbook benchmarks + numerical propagation."""
import math
import sys

from odst import (MU, RE, Mission, bielliptic, evaluate, hohmann, one_tangent,
                  orbit_elements, print_table, simulate, vcirc)

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

print("\n=== 4. Mission case studies ===")
for m in cases:
    print_table(m, evaluate(m))

print(f"\n{'ALL CHECKS PASSED' if fails == 0 else f'{fails} CHECK(S) FAILED'}")
sys.exit(fails)
