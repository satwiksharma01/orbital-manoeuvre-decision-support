"""
Orbital Manoeuvre Decision Support Tool (ODST) - Python reference implementation.

Mirrors the MATLAB implementation in ../matlab exactly. It is used to
  * validate the analytical models against textbook values,
  * verify every manoeuvre's burn plan by numerical propagation (two-body, impulsive),
  * generate the figures used in the report.

Units: km, s, km/s, deg (inputs), rad (internal).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

MU = 398600.4418      # km^3/s^2, Earth
RE = 6378.137         # km
G0 = 9.80665e-3       # km/s^2

# Qualitative operational complexity index (1 = simplest, 5 = hardest).
# Justification is given in docs/REPORT.md section 4.4.
COMPLEXITY = {
    "hohmann": 1.0,
    "hohmann_split": 2.0,
    "one_tangent": 2.5,
    "bielliptic": 3.0,
    "low_thrust": 4.0,
}

WEIGHT_PRESETS = {
    #                    dv     prop   tof   cplx
    "Balanced":          (0.30, 0.20, 0.30, 0.20),
    "Propellant-critical": (0.25, 0.55, 0.05, 0.15),
    "Time-critical":     (0.15, 0.10, 0.60, 0.15),
    "Low-complexity":    (0.20, 0.20, 0.10, 0.50),
}


@dataclass
class Mission:
    alt1: float = 400.0          # initial circular altitude, km
    alt2: float = 35786.0        # target circular altitude, km
    di_deg: float = 28.5         # inclination change, deg
    dv_budget: float = 6.0       # km/s
    tof_max_h: float = 24.0 * 200  # hours
    propulsion: str = "both"     # 'chemical' | 'electric' | 'both'
    isp_chem: float = 320.0      # s
    isp_ep: float = 1600.0       # s
    accel_ep: float = 5e-4       # m/s^2 low-thrust acceleration (F/m)
    rb_max: float = 384400.0     # bi-elliptic apoapsis upper bound, km (lunar distance)
    k_one_tangent: float = 1.5   # one-tangent transfer apsis factor (default)
    weights: tuple = WEIGHT_PRESETS["Balanced"]

    @property
    def r1(self):
        return RE + self.alt1

    @property
    def r2(self):
        return RE + self.alt2


@dataclass
class Option:
    key: str
    name: str
    dv: float                 # total dv, km/s
    tof: float                # s
    burns: int
    isp: float
    complexity: float
    applicable: bool = True
    note: str = ""
    burn_plan: list = field(default_factory=list)  # [(t, v_after, d_inc_rad)]
    params: dict = field(default_factory=dict)
    # filled by the decision layer
    prop_frac: float = float("nan")
    feasible: bool = True
    reasons: list = field(default_factory=list)
    score_ws: float = float("nan")
    score_topsis: float = float("nan")


# ---------------------------------------------------------------- helpers
def vcirc(r):
    return math.sqrt(MU / r)


def vvis(r, a):
    return math.sqrt(MU * (2.0 / r - 1.0 / a))


def combined_burn(va, vb, angle):
    """|dv| between two velocity vectors of magnitude va, vb separated by angle."""
    return math.sqrt(max(va * va + vb * vb - 2.0 * va * vb * math.cos(angle), 0.0))


# ---------------------------------------------------------------- manoeuvres
def hohmann(r1, r2, di, split=None):
    """Hohmann transfer with plane change di split between the two burns.

    split = fraction of di performed at the r1 burn. None -> classical rule:
    whole plane change at the burn on the larger radius (lowest speed).
    """
    if split is None:
        split = 0.0 if r2 >= r1 else 1.0
    if abs(r2 - r1) < 1e-6:            # degenerate: pure plane change at a node
        dv = 2.0 * vcirc(r1) * math.sin(di / 2.0)
        return dv, 0.0, 1, [(0.0, vcirc(r1), di)]
    a = 0.5 * (r1 + r2)
    v1, v2 = vcirc(r1), vcirc(r2)
    vt1, vt2 = vvis(r1, a), vvis(r2, a)
    dv1 = combined_burn(v1, vt1, split * di)
    dv2 = combined_burn(vt2, v2, (1.0 - split) * di)
    tof = math.pi * math.sqrt(a ** 3 / MU)
    plan = [(0.0, vt1, split * di), (tof, v2, (1.0 - split) * di)]
    return dv1 + dv2, tof, 2, plan


def hohmann_optimal_split(r1, r2, di):
    """Numerical optimisation: plane-change split s* minimising total dv."""
    if di <= 0 or abs(r2 - r1) < 1e-6:
        return 0.0 if r2 >= r1 else 1.0
    res = minimize_scalar(lambda s: hohmann(r1, r2, di, s)[0],
                          bounds=(0.0, 1.0), method="bounded",
                          options={"xatol": 1e-8})
    return float(res.x)


def bielliptic(r1, r2, rb, di):
    """Bi-elliptic transfer via apoapsis rb; whole plane change at rb."""
    a1, a2 = 0.5 * (r1 + rb), 0.5 * (r2 + rb)
    dv1 = abs(vvis(r1, a1) - vcirc(r1))
    dv2 = combined_burn(vvis(rb, a1), vvis(rb, a2), di)
    dv3 = abs(vvis(r2, a2) - vcirc(r2))
    t1 = math.pi * math.sqrt(a1 ** 3 / MU)
    t2 = math.pi * math.sqrt(a2 ** 3 / MU)
    plan = [(0.0, vvis(r1, a1), 0.0), (t1, vvis(rb, a2), di), (t1 + t2, vcirc(r2), 0.0)]
    return dv1 + dv2 + dv3, t1 + t2, 3, plan


def bielliptic_optimal_rb(r1, r2, di, rb_max, tof_max=math.inf):
    """Grid search + bounded refinement of rb (dv is multi-modal in rb).
    TOF constraint enforced as a hard bound on rb (TOF increases with rb)."""
    rlo = max(r1, r2) * (1 + 1e-9)
    rhi = rb_max
    if bielliptic(r1, r2, rhi, di)[1] > tof_max:
        if bielliptic(r1, r2, rlo, di)[1] > tof_max:
            return None
        rhi = brentq(lambda rb: bielliptic(r1, r2, rb, di)[1] - tof_max, rlo, rhi)
    grid = np.geomspace(rlo, rhi, 400)
    dvs = np.array([bielliptic(r1, r2, rb, di)[0] for rb in grid])
    i = int(np.argmin(dvs))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    if hi > lo:
        res = minimize_scalar(lambda rb: bielliptic(r1, r2, rb, di)[0],
                              bounds=(lo, hi), method="bounded")
        return float(res.x) if res.fun <= dvs[i] else float(grid[i])
    return float(grid[i])


def one_tangent(r1, r2, k, di):
    """One-tangent-burn (fast) transfer.

    Raising (r2>r1): departure tangential at periapsis r1, transfer apoapsis k*r2 (k>=1).
    Lowering (r2<r1): departure tangential at apoapsis r1, transfer periapsis r2/k.
    The plane change is combined with the burn at the larger radius.
    """
    raising = r2 > r1
    if raising:
        rp, ra = r1, k * r2
    else:
        ra, rp = r1, r2 / k
    a = 0.5 * (rp + ra)
    e = (ra - rp) / (ra + rp)
    p = a * (1 - e * e)
    v1, v2 = vcirc(r1), vcirc(r2)
    vt1, vt2 = vvis(r1, a), vvis(r2, a)
    if e < 1e-12:  # k == 1 and r1 == r2 corner; treat as Hohmann
        return hohmann(r1, r2, di)
    cosnu = max(-1.0, min(1.0, (p / r2 - 1.0) / e))
    nu = math.acos(cosnu)                         # 0..pi, symmetric arc
    phi = math.atan2(e * math.sin(nu), 1 + e * math.cos(nu))  # flight-path angle
    E = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))
    t_from_peri = math.sqrt(a ** 3 / MU) * (E - e * math.sin(E))
    if raising:
        tof = t_from_peri
        dv1 = abs(vt1 - v1)
        dv2 = math.sqrt(max(vt2 ** 2 + v2 ** 2 - 2 * vt2 * v2 * math.cos(phi) * math.cos(di), 0))
        plan = [(0.0, vt1, 0.0), (tof, v2, di)]
    else:
        tof = math.pi * math.sqrt(a ** 3 / MU) - t_from_peri
        dv1 = combined_burn(v1, vt1, di)
        dv2 = math.sqrt(max(vt2 ** 2 + v2 ** 2 - 2 * vt2 * v2 * math.cos(phi), 0))
        plan = [(0.0, vt1, di), (tof, v2, 0.0)]
    return dv1 + dv2, tof, 2, plan, {"e": e, "nu_deg": math.degrees(nu), "fpa_deg": math.degrees(phi)}


def one_tangent_for_tof(r1, r2, di, tof_target, kmax=50.0):
    """Smallest k (=> minimum dv) whose one-tangent TOF meets tof_target."""
    f = lambda k: one_tangent(r1, r2, k, di)[1] - tof_target
    if f(1.0 + 1e-9) <= 0:
        return 1.0 + 1e-9
    if f(kmax) > 0:
        return None
    return brentq(f, 1.0 + 1e-9, kmax, xtol=1e-10)


def low_thrust(r1, r2, di, accel_ms2):
    """Edelbaum continuous low-thrust transfer between circular orbits."""
    v1, v2 = vcirc(r1), vcirc(r2)
    dv = math.sqrt(max(v1 ** 2 + v2 ** 2 - 2 * v1 * v2 * math.cos(math.pi / 2 * di), 0))
    tof = dv / (accel_ms2 * 1e-3)
    return dv, tof


# ---------------------------------------------------------------- evaluation
def evaluate(m: Mission):
    r1, r2, di = m.r1, m.r2, math.radians(m.di_deg)
    tof_max = m.tof_max_h * 3600.0
    opts = []

    dv, tof, nb, plan = hohmann(r1, r2, di)
    coplanar_same = abs(r2 - r1) < 1e-6
    opts.append(Option("hohmann", "Simple plane change" if coplanar_same else "Hohmann (Δi at high burn)",
                       dv, tof, nb, m.isp_chem, COMPLEXITY["hohmann"], burn_plan=plan))

    s = hohmann_optimal_split(r1, r2, di)
    dv, tof, nb, plan = hohmann(r1, r2, di, s)
    o = Option("hohmann_split", "Hohmann (optimal Δi split)", dv, tof, nb, m.isp_chem,
               COMPLEXITY["hohmann_split"], burn_plan=plan, params={"split_at_r1": s})
    if di <= 0 or coplanar_same:
        o.applicable, o.note = False, "identical to Hohmann (no plane change to split)"
    opts.append(o)

    rb = bielliptic_optimal_rb(r1, r2, di, max(m.rb_max, 1.01 * max(r1, r2)), tof_max)
    tof_ok = rb is not None
    if not tof_ok:                      # even the shortest bi-elliptic misses the deadline
        rb = 2.0 * max(r1, r2)
    dv, tof, nb, plan = bielliptic(r1, r2, rb, di)
    o = Option("bielliptic", "Bi-elliptic", dv, tof, nb, m.isp_chem, COMPLEXITY["bielliptic"],
               burn_plan=plan, params={"rb": rb, "rb_over_r1": rb / r1})
    if tof_ok and rb < max(r1, r2) * 1.001:
        o.applicable, o.note = False, "optimum degenerates to Hohmann (rB -> r2)"
    opts.append(o)

    if not coplanar_same:
        tof_h = hohmann(r1, r2, di)[1]
        if tof_max < tof_h:
            k = one_tangent_for_tof(r1, r2, di, tof_max)
            note = "k solved to meet TOF limit" if k else "cannot meet TOF limit (k>50)"
            k = k or 50.0
        else:
            k, note = m.k_one_tangent, "default k"
        dv, tof, nb, plan, extra = one_tangent(r1, r2, k, di)
        extra.update({"k": k})
        opts.append(Option("one_tangent", "One-tangent burn (fast)", dv, tof, nb, m.isp_chem,
                           COMPLEXITY["one_tangent"], burn_plan=plan, params=extra, note=note))

    dv, tof = low_thrust(r1, r2, di, m.accel_ep)
    opts.append(Option("low_thrust", "Low-thrust spiral (Edelbaum)", dv, tof, 0, m.isp_ep,
                       COMPLEXITY["low_thrust"]))

    for o in opts:
        o.prop_frac = 1.0 - math.exp(-o.dv / (o.isp * G0))
    apply_rules(m, opts)
    score(m, opts)
    return opts



# ------------------------------------------------- rule-based decision layer
def apply_rules(m: Mission, opts):
    tof_max = m.tof_max_h * 3600.0
    for o in opts:
        o.feasible, o.reasons = o.applicable, ([] if o.applicable else [o.note])
        if not o.applicable:
            continue
        if o.dv > m.dv_budget:
            o.feasible = False
            o.reasons.append(f"Δv {o.dv:.3f} > budget {m.dv_budget:.3f} km/s")
        if o.tof > tof_max * (1 + 1e-6):
            o.feasible = False
            o.reasons.append(f"TOF {o.tof/3600:.1f} h > limit {m.tof_max_h:.1f} h")
        if o.key == "low_thrust" and m.propulsion == "chemical":
            o.feasible = False
            o.reasons.append("no electric propulsion on board")
        if o.key != "low_thrust" and m.propulsion == "electric":
            o.feasible = False
            o.reasons.append("no high-thrust (chemical) propulsion on board")


def rule_advice(m: Mission):
    """Textbook rules of thumb used as the 'rule-based decision matrix'."""
    r1, r2 = m.r1, m.r2
    R = max(r1, r2) / min(r1, r2)
    adv = []
    if R < 11.94:
        adv.append(f"R = {R:.2f} < 11.94: Hohmann is the minimum-Δv two/three-impulse transfer.")
    elif R < 15.58:
        adv.append(f"11.94 ≤ R = {R:.2f} < 15.58: bi-elliptic wins only if rB is large enough.")
    else:
        adv.append(f"R = {R:.2f} ≥ 15.58: bi-elliptic beats Hohmann for any rB > r2.")
    if m.di_deg > 0:
        adv.append("Plane change: perform it where speed is lowest (largest radius) or combine with a burn.")
    if R < 1.0 + 1e-9 and m.di_deg > 48.94:
        adv.append("Δi > 48.94°: a three-impulse (bi-elliptic) plane change beats a single burn.")
    tof_h = hohmann(r1, r2, math.radians(m.di_deg))[1] / 3600
    if m.tof_max_h < tof_h:
        adv.append(f"TOF limit {m.tof_max_h:.2f} h < Hohmann {tof_h:.2f} h: a fast (one-tangent) transfer is required.")
    if m.propulsion != "chemical":
        adv.append("Electric propulsion available: low-thrust spiral minimises propellant if time allows.")
    return adv


# ------------------------------------------------------------ MCDA scoring
def _criteria(o):
    # TOF on a log scale: hours vs months differ in kind, not linearly.
    return np.array([o.dv, o.prop_frac, math.log10(max(o.tof, 60.0)), o.complexity])


def score(m: Mission, opts):
    feas = [o for o in opts if o.feasible]
    if not feas:
        return
    w = np.asarray(m.weights, float)
    w = w / w.sum()
    X = np.array([_criteria(o) for o in feas])
    # (1) weighted sum, min-max normalised, all criteria are costs -> 1 = best
    lo, hi = X.min(0), X.max(0)
    span = np.where(hi - lo > 1e-12, hi - lo, 1.0)
    U = np.where(hi - lo > 1e-12, (hi - X) / span, 1.0)
    ws = U @ w
    # (2) TOPSIS, vector normalisation
    nrm = np.linalg.norm(X, axis=0)
    nrm[nrm < 1e-12] = 1.0
    V = X / nrm * w
    best, worst = V.min(0), V.max(0)
    dp = np.linalg.norm(V - best, axis=1)
    dm = np.linalg.norm(V - worst, axis=1)
    cc = np.where(dp + dm > 1e-12, dm / np.maximum(dp + dm, 1e-12), 1.0)
    for o, a, b in zip(feas, ws, cc):
        o.score_ws, o.score_topsis = float(a), float(b)


def recommend(opts):
    feas = [o for o in opts if o.feasible]
    if not feas:
        return None
    return max(feas, key=lambda o: o.score_ws)


# ------------------------------------------------------ numerical verification
def _twobody(t, y):
    r = y[:3]
    return np.concatenate([y[3:], -MU * r / np.linalg.norm(r) ** 3])


def _rotate(v, axis, ang):
    """Rodrigues rotation of v about unit axis by ang."""
    k = axis / np.linalg.norm(axis)
    return v * math.cos(ang) + np.cross(k, v) * math.sin(ang) + k * np.dot(k, v) * (1 - math.cos(ang))


def simulate(r1, plan, coast_after=0.0, n_out=2000):
    """Propagate an impulsive burn plan from a circular equatorial orbit at r1.
    Every burn in the ODST plans leaves velocity perpendicular to r (apsis or
    circularisation), so a burn = (time, speed after, plane rotation about r_hat)."""
    y = np.array([r1, 0, 0, 0, vcirc(r1), 0.0])
    t0, ts, ys, dv_used = 0.0, [], [], 0.0
    times = [b[0] for b in plan] + [plan[-1][0] + coast_after]
    for i, (tb, vafter, dinc) in enumerate(plan):
        if tb > t0:
            sol = solve_ivp(_twobody, (t0, tb), y, rtol=1e-11, atol=1e-11, dense_output=True)
            tt = np.linspace(t0, tb, max(int(n_out * (tb - t0) / times[-1]), 50))
            ts.append(tt); ys.append(sol.sol(tt).T)
            y = sol.y[:, -1].copy()
        r, v = y[:3], y[3:]
        h = np.cross(r, v)
        # rotate about the node line in the sense that increases inclination
        # w.r.t. the initial plane (keeps split plane changes additive)
        nh = h / np.linalg.norm(h)
        cands = [_rotate(nh, r, dinc), _rotate(nh, r, -dinc)]
        n_new = min(cands, key=lambda n: n[2])
        vnew = vafter * np.cross(n_new, r / np.linalg.norm(r))
        dv_used += np.linalg.norm(vnew - v)
        y = np.concatenate([r, vnew])
        t0 = tb
    if coast_after > 0:
        sol = solve_ivp(_twobody, (t0, t0 + coast_after), y, rtol=1e-11, atol=1e-11, dense_output=True)
        tt = np.linspace(t0, t0 + coast_after, 400)
        ts.append(tt); ys.append(sol.sol(tt).T)
        y = sol.y[:, -1]
    return np.concatenate(ts), np.vstack(ys), y, dv_used


def orbit_elements(y):
    r, v = y[:3], y[3:]
    h = np.cross(r, v)
    ev = np.cross(v, h) / MU - r / np.linalg.norm(r)
    a = 1.0 / (2 / np.linalg.norm(r) - np.dot(v, v) / MU)
    inc = math.degrees(math.acos(h[2] / np.linalg.norm(h)))
    return a, np.linalg.norm(ev), inc


# ------------------------------------------------------------------- report
def print_table(m: Mission, opts):
    print(f"\nMission: alt {m.alt1:.0f} -> {m.alt2:.0f} km, Δi = {m.di_deg:.2f}°, "
          f"budget {m.dv_budget} km/s, TOF max {m.tof_max_h:.1f} h, weights {m.weights}")
    hdr = f"{'Option':34s}{'Δv km/s':>9s}{'TOF':>12s}{'burns':>6s}{'mp/m0':>8s}{'cplx':>6s}{'WS':>7s}{'TOPSIS':>8s}  status"
    print(hdr); print("-" * len(hdr))
    for o in opts:
        tof = f"{o.tof/3600:.2f} h" if o.tof < 72 * 3600 else f"{o.tof/86400:.1f} d"
        st = "OK" if o.feasible else "X " + "; ".join(o.reasons)
        ws = f"{o.score_ws:.3f}" if o.feasible else "-"
        tp = f"{o.score_topsis:.3f}" if o.feasible else "-"
        print(f"{o.name:34s}{o.dv:9.4f}{tof:>12s}{o.burns:6d}{o.prop_frac:8.3f}{o.complexity:6.1f}{ws:>7s}{tp:>8s}  {st}")
    best = recommend(opts)
    print("Recommended:", best.name if best else "NONE FEASIBLE - relax constraints")
    for a in rule_advice(m):
        print("  rule:", a)


if __name__ == "__main__":
    print_table(Mission(), evaluate(Mission()))
