"""
Orbital Manoeuvre Decision Support Tool (ODST) - Python reference implementation.

Mirrors the MATLAB implementation in ../matlab exactly. It is used to
  * validate the analytical models against textbook values,
  * verify every manoeuvre's burn plan by numerical propagation (two-body, impulsive),
  * generate the figures used in the report,
  * export a transfer's burn plan as a NASA GMAT script and as JSON.

Units: km, s, km/s, deg (inputs), rad (internal).
"""
from __future__ import annotations

import json
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
    "apse": 1.5,
    "lambert": 2.5,
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
    alt1: float = 400.0          # initial altitude (periapsis if elliptical), km
    alt2: float = 35786.0        # target altitude (periapsis if elliptical), km
    apo1: float | None = None    # initial apoapsis altitude, km; None = circular
    apo2: float | None = None    # target apoapsis altitude, km; None = circular
    dw_deg: float = 0.0          # target periapsis, deg from the node line (the start periapsis if w1_deg = 0)
    di_deg: float = 28.5         # inclination change (angle between the orbit planes), deg
    dv_budget: float = 6.5       # km/s
    margin_pct: float = 5.0      # Δv margin the budget must also cover, %
    tof_max_h: float = 24.0 * 200  # hours
    propulsion: str = "both"     # 'chemical' | 'electric' | 'both'
    isp_chem: float = 320.0      # s
    tw_chem: float = 0.5         # chemical thrust-to-weight at departure, T / (m0 g0)
    isp_ep: float = 1600.0       # s
    accel_ep: float = 5e-4       # m/s^2 low-thrust acceleration (F/m)
    rb_max: float = 384400.0     # bi-elliptic apoapsis upper bound, km (lunar distance)
    k_one_tangent: float = 1.5   # one-tangent transfer apsis factor (default)
    weights: tuple = WEIGHT_PRESETS["Balanced"]
    w1_deg: float = 0.0          # start periapsis, deg from the node line
    elements: tuple | None = None  # (i1, RAAN1, argp1, i2, RAAN2, argp2), deg; sets di_deg, w1_deg, dw_deg
    R: object = field(default=None, init=False, repr=False, compare=False)  # relative frame -> EarthMJ2000Eq
    # spacecraft model: with mass0_kg set, the masses, thrust and power below set tw_chem and accel_ep,
    # the electric thrust is constant (the spiral speeds up as the mass falls) and the propellant on board is a limit
    mass0_kg: float | None = None  # wet mass at departure, kg; None = the generic T/W and F/m inputs
    thrust_chem_N: float = 4900.0  # chemical engine thrust, N
    prop_chem_kg: float = 800.0    # chemical propellant on board, kg
    power_ep_kW: float = 6.5       # electric thruster input power, kW
    eff_ep: float = 0.6            # electric thruster efficiency: jet power / input power
    prop_ep_kg: float = 400.0      # electric propellant on board, kg

    def __post_init__(self):
        if self.elements is not None:
            self.di_deg, self.w1_deg, self.dw_deg, self.R = relative_geometry(self)
        if self.mass0_kg is not None:
            self.tw_chem = self.thrust_chem_N / (self.mass0_kg * G0 * 1e3)
            self.accel_ep = ep_thrust(self) / self.mass0_kg

    @property
    def r1(self):
        return RE + self.alt1

    @property
    def r2(self):
        return RE + self.alt2

    @property
    def ra1(self):
        return RE + (self.alt1 if self.apo1 is None else self.apo1)

    @property
    def ra2(self):
        return RE + (self.alt2 if self.apo2 is None else self.apo2)

    @property
    def elliptical(self):
        return self.ra1 - self.r1 > ELL_TOL or self.ra2 - self.r2 > ELL_TOL


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
    burn_plan: list = field(default_factory=list)  # [(t, v_after, d_inc_rad, r_burn, dv_burn)]
    vec_plan: list = field(default_factory=list)   # every impulsive option: [(t, dv_vec)] from state y0
    y0: object = None
    legs: list = field(default_factory=list)       # [(r_burn, dv_burn)] for the gravity loss
    params: dict = field(default_factory=dict)
    # filled by the decision layer
    dv_loss: float = 0.0      # finite-burn gravity loss, km/s
    prop_frac: float = float("nan")
    prop_kg: float = float("nan")  # propellant mass, kg, with the spacecraft model
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
        return dv, 0.0, 1, [(0.0, vcirc(r1), di, r1, dv)]
    a = 0.5 * (r1 + r2)
    v1, v2 = vcirc(r1), vcirc(r2)
    vt1, vt2 = vvis(r1, a), vvis(r2, a)
    dv1 = combined_burn(v1, vt1, split * di)
    dv2 = combined_burn(vt2, v2, (1.0 - split) * di)
    tof = math.pi * math.sqrt(a ** 3 / MU)
    plan = [(0.0, vt1, split * di, r1, dv1), (tof, v2, (1.0 - split) * di, r2, dv2)]
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
    plan = [(0.0, vvis(r1, a1), 0.0, r1, dv1), (t1, vvis(rb, a2), di, rb, dv2), (t1 + t2, vcirc(r2), 0.0, r2, dv3)]
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
        plan = [(0.0, vt1, 0.0, r1, dv1), (tof, v2, di, r2, dv2)]
    else:
        tof = math.pi * math.sqrt(a ** 3 / MU) - t_from_peri
        dv1 = combined_burn(v1, vt1, di)
        dv2 = math.sqrt(max(vt2 ** 2 + v2 ** 2 - 2 * vt2 * v2 * math.cos(phi), 0))
        plan = [(0.0, vt1, di, r1, dv1), (tof, v2, 0.0, r2, dv2)]
    return dv1 + dv2, tof, 2, plan, {"e": e, "nu_deg": math.degrees(nu), "fpa_deg": math.degrees(phi)}


def one_tangent_for_tof(r1, r2, di, tof_target, kmax=50.0):
    """Smallest k (=> minimum dv) whose one-tangent TOF meets tof_target."""
    f = lambda k: one_tangent(r1, r2, k, di)[1] - tof_target
    if f(1.0 + 1e-9) <= 0:
        return 1.0 + 1e-9
    if f(kmax) > 0:
        return None
    return brentq(f, 1.0 + 1e-9, kmax, xtol=1e-10)


EDELBAUM_MAX_DI = 2.0   # rad (114.59 deg): the yaw angle sweeps pi/2 * di, which must not exceed pi


def low_thrust(r1, r2, di, accel_ms2, isp=None):
    """Edelbaum continuous low-thrust transfer between circular orbits, valid for di <= 2 rad.
    Beyond that the closed form's dv falls as di grows and its yaw law misses the target plane.
    accel_ms2 is the acceleration at departure. Without isp it stays constant: tof = dv / a. With
    isp the thrust is constant and the mass falls: tof = ve / a (1 - exp(-dv / ve)). Edelbaum's dv
    and yaw law are unchanged, as both depend only on the velocity increment u = integral of F/m dt."""
    v1, v2 = vcirc(r1), vcirc(r2)
    dv = math.sqrt(max(v1 ** 2 + v2 ** 2 - 2 * v1 * v2 * math.cos(math.pi / 2 * di), 0))
    if isp is None:
        return dv, dv / (accel_ms2 * 1e-3)
    ve = isp * G0
    return dv, ve / (accel_ms2 * 1e-3) * (1 - math.exp(-dv / ve))


def ep_thrust(m):
    """Electric thrust, N, from input power: F = 2 eta P / ve (the jet power is eta P)."""
    return 2 * m.eff_ep * m.power_ep_kW * 1e3 / (m.isp_ep * G0 * 1e3)


# ------------------------------------------- elliptical orbits and Lambert
# Relative frame: the start orbit lies in the x-y plane and +x is the node line, where the
# target crosses that plane going north (+z). The target plane is turned by di about +x.
# Each periapsis is measured from +x in its own plane: w1 for the start orbit (0 unless
# set from orbital elements), dw for the target. Burns are vectors here, so these plans
# are (t, dv_vec) rows with the start state y0, not the circular format.
ELL_TOL = 1e-6          # km: apoapsis - periapsis below this is a circular orbit
DV_ZERO = 1e-9          # km/s: burns smaller than this are dropped
H_ARC_MIN = 100.0       # km: a Lambert transfer arc may not pass below this altitude
K_SEARCH, K_POLISH = 0.01, 10.0   # km/s per km below H_ARC_MIN: gentle penalty to search, stiff to finish


def _on_node_line(ang):
    return abs(math.sin(ang)) < 1e-9


def orbit_basis(di, dw):
    """Unit vectors to periapsis (P) and 90 deg ahead of it (Q) for a plane turned by
    di about +x with periapsis dw from +x."""
    x, y = np.array([1.0, 0.0, 0.0]), np.array([0.0, math.cos(di), math.sin(di)])
    return math.cos(dw) * x + math.sin(dw) * y, -math.sin(dw) * x + math.cos(dw) * y


def perifocal(inc, raan, argp):
    """Unit vectors to periapsis (P), 90 deg ahead of it (Q) and along the orbit normal (W)
    in EarthMJ2000Eq, from inclination, RAAN and argument of periapsis (rad)."""
    ci, si, cO, sO, cw, sw = math.cos(inc), math.sin(inc), math.cos(raan), math.sin(raan), math.cos(argp), math.sin(argp)
    return (np.array([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si]),
            np.array([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si]),
            np.array([sO * si, -cO * si, ci]))


def relative_geometry(m):
    """(di_deg, w1_deg, dw_deg, R) for a mission given by orbital elements: the angle between
    the planes, each periapsis from the node line, and R, whose columns are the relative
    frame's axes in EarthMJ2000Eq. Coplanar orbits have no node line, so the start periapsis
    takes its place; a circular orbit's periapsis is put on the node line."""
    P1, _, W1 = perifocal(*map(math.radians, m.elements[:3]))
    P2, _, W2 = perifocal(*map(math.radians, m.elements[3:]))
    n = np.cross(W1, W2)
    x = n / np.linalg.norm(n) if np.linalg.norm(n) > 1e-9 else P1
    y, y2 = np.cross(W1, x), np.cross(W2, x)                  # 90 deg ahead of +x in each plane
    ang = lambda P, yy, circ: 0.0 if circ else round(math.degrees(math.atan2(P @ yy, P @ x)), 12) % 360   # -1e-16 -> 0, not 360
    return (math.degrees(math.atan2(np.linalg.norm(n), W1 @ W2)), ang(P1, y, m.ra1 - m.r1 < ELL_TOL),
            ang(P2, y2, m.ra2 - m.r2 < ELL_TOL), np.column_stack([x, y, W1]))


def orbit_radius(rp, ra, th):
    if abs(math.cos(th) - 1) < 1e-12:
        return rp
    if abs(math.cos(th) + 1) < 1e-12:
        return ra
    e, p = (ra - rp) / (ra + rp), 2 * rp * ra / (rp + ra)
    return p / (1 + e * math.cos(th))


def orbit_state(rp, ra, P, Q, th):
    """Position and velocity at true anomaly th."""
    e, p = (ra - rp) / (ra + rp), 2 * rp * ra / (rp + ra)
    r = orbit_radius(rp, ra, th)
    return (r * (math.cos(th) * P + math.sin(th) * Q),
            math.sqrt(MU / p) * (-math.sin(th) * P + (e + math.cos(th)) * Q))


def arc_min_radius(r1v, v1v, r2v):
    """Lowest radius on the conic from r1v (velocity v1v) forward to r2v: the transfer's
    periapsis radius if periapsis lies between them, otherwise the lower end point."""
    h = np.cross(r1v, v1v)
    ev = np.cross(v1v, h) / MU - r1v / np.linalg.norm(r1v)
    e, ends = np.linalg.norm(ev), min(np.linalg.norm(r1v), np.linalg.norm(r2v))
    if e < 1e-10:
        return ends
    q = np.cross(h, ev)                      # 90 deg ahead of periapsis, in the direction of motion
    nu = lambda r: math.atan2(np.dot(q, r) / np.linalg.norm(q), np.dot(ev, r) / e) % (2 * math.pi)
    if nu(r2v) < nu(r1v):                    # moving forward from r1 to r2 wraps through periapsis
        return min(ends, np.dot(h, h) / MU / (1 + e))
    return ends


def stumpff(z):
    if z > 1e-8:
        s = math.sqrt(z)
        return (1 - math.cos(s)) / z, (s - math.sin(s)) / s ** 3
    if z < -1e-8:
        s = math.sqrt(-z)
        return (math.cosh(s) - 1) / -z, (math.sinh(s) - s) / s ** 3
    return 0.5 - z / 24, 1 / 6 - z / 120


def lambert(r1v, r2v, dt):
    """Single-revolution prograde Lambert problem (transfer angular momentum along +z),
    universal variables (Curtis Algorithm 5.2) with bisection on z (Vallado's
    bracketing). Returns (v1, v2), or None if there is no solution."""
    r1, r2 = np.linalg.norm(r1v), np.linalg.norm(r2v)
    th = math.acos(max(-1.0, min(1.0, np.dot(r1v, r2v) / (r1 * r2))))
    if np.cross(r1v, r2v)[2] < 0:
        th = 2 * math.pi - th
    if 1 - math.cos(th) < 1e-14:
        return None
    A = math.sin(th) * math.sqrt(r1 * r2 / (1 - math.cos(th)))
    if abs(A) < 1e-12:
        return None

    def y_of(z):
        C, S = stumpff(z)
        return r1 + r2 + A * (z * S - 1) / math.sqrt(C), C, S

    def tof(z):
        y, C, S = y_of(z)
        if y < 0:
            return -math.inf
        x = math.sqrt(y / C)
        return (x ** 3 * S + A * math.sqrt(y)) / math.sqrt(MU)

    lo, hi = -4 * math.pi, 4 * math.pi ** 2 - 1e-4   # C(z) -> 0 at 4 pi^2
    while tof(lo) > dt:                       # very fast: needs a more hyperbolic bracket
        lo *= 2
        if lo < -1e5:
            return None
    if tof(hi) < dt:
        return None
    for _ in range(200):
        z = 0.5 * (lo + hi)
        if tof(z) < dt:
            lo = z
        else:
            hi = z
        if hi - lo < 1e-13 * max(1.0, abs(z)):
            break
    y = y_of(0.5 * (lo + hi))[0]
    f, g, gdot = 1 - y / r1, A * math.sqrt(y / MU), 1 - y / r2
    return (r2v - f * r1v) / g, (gdot * r2v - r1v) / g


def apse_transfer(rp1, ra1, rp2, ra2, di, dw, w1=0.0):
    """Two-impulse transfer from an apsis of the start orbit to the opposite direction
    on the target orbit, with the plane change split optimally between the burns.
    Returns (dict, "") or (None, reason)."""
    circ1, circ2 = ra1 - rp1 < ELL_TOL, ra2 - rp2 < ELL_TOL
    w1 = 0.0 if circ1 else w1                     # a circular orbit has no periapsis to place
    if not (circ1 or circ2 or _on_node_line(dw - w1)):
        return None, f"apse lines {math.degrees((dw - w1) % (2 * math.pi)):g}° apart, not 0° or 180°"
    dirs = [dw, dw + math.pi] if circ1 and not circ2 else [w1, w1 + math.pi]   # departure, from +x
    if di > 0 and not all(_on_node_line(d) for d in dirs):
        return None, "a plane change needs both burns on the node line, and the apse line isn't on it"
    P1, Q1 = orbit_basis(0.0, w1)
    P2, Q2 = orbit_basis(di, dw)
    a1, a2 = (rp1 + ra1) / 2, (rp2 + ra2) / 2
    best = None
    for phi in dirs:
        th2 = phi + math.pi - dw
        rA, rB = orbit_radius(rp1, ra1, phi - w1), orbit_radius(rp2, ra2, th2)
        at = (rA + rB) / 2
        v1, vtA, vtB, v2 = vvis(rA, a1), vvis(rA, at), vvis(rB, at), vvis(rB, a2)
        cost = lambda s: combined_burn(v1, vtA, s * di) + combined_burn(vtB, v2, (1 - s) * di)
        cands = [0.0, 1.0]
        if di > 0:
            cands.append(float(minimize_scalar(cost, bounds=(0, 1), method="bounded",
                                               options={"xatol": 1e-10}).x))
        s = min(cands, key=cost)
        if best is None or cost(s) < best[0] - 1e-12:
            best = (cost(s), phi, th2, rA, rB, at, v1, vtA, vtB, v2, s)
    _, phi, th2, rA, rB, at, v1, vtA, vtB, v2, s = best
    dv1, dv2 = combined_burn(v1, vtA, s * di), combined_burn(vtB, v2, (1 - s) * di)
    Rt = np.array([[1, 0, 0], [0, math.cos(s * di), -math.sin(s * di)], [0, math.sin(s * di), math.cos(s * di)]])
    along = lambda ang: Rt @ np.array([-math.sin(ang), math.cos(ang), 0.0])
    r1v, v1v = orbit_state(rp1, ra1, P1, Q1, phi - w1)
    rBv = Rt @ np.array([rB * math.cos(phi + math.pi), rB * math.sin(phi + math.pi), 0.0])
    vB = vtB * along(phi + math.pi)
    dv1v = vtA * along(phi) - v1v
    dv2v = orbit_state(rp2, ra2, P2, Q2, th2)[1] - vB
    tof = math.pi * math.sqrt(at ** 3 / MU)
    keep1, keep2 = dv1 > DV_ZERO, dv2 > DV_ZERO
    if keep1 and keep2:
        y0, plan, legs = np.concatenate([r1v, v1v]), [(0.0, dv1v), (tof, dv2v)], [(rA, dv1), (rB, dv2)]
    elif keep2:                                   # already on the transfer orbit: one burn on arrival
        y0, plan, legs, tof = np.concatenate([rBv, vB]), [(0.0, dv2v)], [(rB, dv2)], 0.0
    else:                                         # the first burn lands on the target orbit
        y0, plan, legs, tof = np.concatenate([r1v, v1v]), [(0.0, dv1v)] if keep1 else [], [(rA, dv1)] if keep1 else [], 0.0
    return {"dv": dv1 * keep1 + dv2 * keep2, "tof": tof, "burns": len(plan), "y0": y0, "plan": plan,
            "legs": legs, "split_at_dep": s, "depart_deg": math.degrees(phi - w1) % 360,
            "rA": rA, "rB": rB}, ""


def _nelder_mead(f, x0, steps, iters=400, tol=1e-12):
    """Plain Nelder-Mead (reflect 1, expand 2, contract 0.5, shrink 0.5), written out so
    the MATLAB and web engines can run the identical search."""
    pts = [list(x0)] + [[x0[j] + (steps[j] if j == i else 0.0) for j in range(len(x0))] for i in range(len(x0))]
    vals = [f(p) for p in pts]
    for _ in range(iters):
        order = sorted(range(len(pts)), key=lambda i: vals[i])
        pts, vals = [pts[i] for i in order], [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol:
            break
        n = len(x0)
        c = [sum(p[j] for p in pts[:-1]) / n for j in range(n)]
        xr = [c[j] + (c[j] - pts[-1][j]) for j in range(n)]
        fr = f(xr)
        if fr < vals[0]:
            xe = [c[j] + 2 * (c[j] - pts[-1][j]) for j in range(n)]
            fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = [c[j] + 0.5 * (pts[-1][j] - c[j]) for j in range(n)]
            fc = f(xc)
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                pts = [pts[0]] + [[pts[0][j] + 0.5 * (p[j] - pts[0][j]) for j in range(n)] for p in pts[1:]]
                vals = [vals[0]] + [f(p) for p in pts[1:]]
    i = min(range(len(pts)), key=lambda k: vals[k])
    return pts[i], vals[i]


def lambert_optimal(rp1, ra1, rp2, ra2, di, dw, tof_max, n_grid=16, w1=0.0):
    """Minimum-dv two-impulse Lambert transfer over departure true anomaly, arrival true
    anomaly and transfer time (60 s to min(tof_max, the longer orbital period)).
    Transfers within 0.5 deg of 0 or 180 deg are skipped: there the transfer plane is
    undefined. Grid n_grid x n_grid with a golden-section search on time per cell, then
    Nelder-Mead from the best cell, restarted once.
    The arc must stay above H_ARC_MIN altitude. That floor is an exact penalty, not a wall:
    a wall stops Nelder-Mead short of the boundary (46 m/s high on a floor-active case in the
    2026-10-22 audit). The search uses K_SEARCH and the final polish K_POLISH, so the result
    sits on or above the floor; a result still below it is rejected."""
    P1, Q1 = orbit_basis(0.0, w1)
    P2, Q2 = orbit_basis(di, dw)
    tmin = 60.0
    tmax = min(tof_max, max(2 * math.pi * math.sqrt(((rp + ra) / 2) ** 3 / MU) for rp, ra in ((rp1, ra1), (rp2, ra2))))
    if tmax <= tmin:
        return None
    gap = math.sin(math.radians(0.5))

    def solve(x, K):
        th1, th2, u = x
        if not 0.0 <= u <= 1.0:
            return math.inf, None
        r1v, v1v = orbit_state(rp1, ra1, P1, Q1, th1)
        r2v, v2v = orbit_state(rp2, ra2, P2, Q2, th2)
        if np.linalg.norm(np.cross(r1v, r2v)) < gap * np.linalg.norm(r1v) * np.linalg.norm(r2v):
            return math.inf, None
        sol = lambert(r1v, r2v, tmin + u * (tmax - tmin))
        if sol is None:
            return math.inf, None
        dv1v, dv2v = sol[0] - v1v, v2v - sol[1]
        dv = np.linalg.norm(dv1v) + np.linalg.norm(dv2v)
        short = RE + H_ARC_MIN - arc_min_radius(r1v, sol[0], r2v)          # km below the floor, if > 0
        return dv + K * max(0.0, short), (r1v, v1v, r2v, dv1v, dv2v, dv, short)

    cost = lambda x: solve(x, K_SEARCH)[0]
    g = (math.sqrt(5) - 1) / 2
    best = (math.inf, None)
    for i in range(n_grid):
        for j in range(n_grid):
            th1, th2 = 2 * math.pi * i / n_grid, 2 * math.pi * j / n_grid
            a, b = 0.0, 1.0
            c, d = b - g * (b - a), a + g * (b - a)
            fc, fd = cost((th1, th2, c)), cost((th1, th2, d))
            for _ in range(25):
                if fc < fd:
                    b, d, fd = d, c, fc
                    c = b - g * (b - a)
                    fc = cost((th1, th2, c))
                else:
                    a, c, fc = c, d, fd
                    d = a + g * (b - a)
                    fd = cost((th1, th2, d))
            u = (a + b) / 2
            fu = cost((th1, th2, u))
            if fu < best[0]:
                best = (fu, [th1, th2, u])
    if not math.isfinite(best[0]):
        return None
    step = 2 * math.pi / n_grid
    x, _ = _nelder_mead(cost, best[1], [step / 2, step / 2, 0.05])
    x, _ = _nelder_mead(cost, x, [step / 4, step / 4, 0.02])                 # restart: Nelder-Mead can stall
    x, _ = _nelder_mead(lambda y: solve(y, K_POLISH)[0], x, [step / 8, step / 8, 0.01])
    _, (r1v, v1v, r2v, dv1v, dv2v, fx, short) = solve(x, K_POLISH)
    if short > 1e-6:
        return None
    T = tmin + x[2] * (tmax - tmin)
    return {"dv": fx, "tof": T, "burns": 2, "y0": np.concatenate([r1v, v1v]), "plan": [(0.0, dv1v), (T, dv2v)],
            "legs": [(np.linalg.norm(r1v), np.linalg.norm(dv1v)), (np.linalg.norm(r2v), np.linalg.norm(dv2v))],
            "th1_deg": math.degrees(x[0]) % 360, "th2_deg": math.degrees(x[1]) % 360,
            "time_limited": x[2] > 1 - 1e-6}


def gravity_loss(legs, isp, tw):
    """Finite-burn gravity loss, Robbins (1966): dv * (w tau)^2 / 24 per burn, with
    w the circular rate at the burn radius and tau the burn time at constant thrust
    T = tw * m0 * g0. Matches a fixed-attitude burn within ~5 % for tw >= 0.2 and
    overestimates a velocity-tracking burn about 2x (validate.py, section 4)."""
    ve, m, loss = isp * G0, 1.0, 0.0
    for r, dv in legs:
        tau = isp / tw * m * (1.0 - math.exp(-dv / ve))
        loss += dv * (math.sqrt(MU / r ** 3) * tau) ** 2 / 24.0
        m *= math.exp(-dv / ve)
    return loss


# ---------------------------------------------------------------- evaluation
def evaluate(m: Mission):
    opts = _elliptical_options(m) if m.elliptical else _circular_options(m)
    for o in opts:
        if o.burn_plan:
            o.legs = [(b[3], b[4]) for b in o.burn_plan]
            o.y0, o.vec_plan = circular_vec_plan(m.r1, o.burn_plan, math.radians(m.di_deg))
        if o.legs:
            o.dv_loss = gravity_loss(o.legs, o.isp, m.tw_chem)
        # the propellant also pays for the finite-burn loss, as the budget check does
        o.prop_frac = 1.0 - math.exp(-(o.dv + o.dv_loss) / (o.isp * G0))
        if m.mass0_kg is not None:
            o.prop_kg = m.mass0_kg * o.prop_frac
    apply_rules(m, opts)
    score(m, opts)
    return opts


def _elliptical_options(m: Mission):
    """Hohmann, bi-elliptic, one-tangent and Edelbaum assume circular orbits, so an
    elliptical mission is sized with the apse-to-apse and Lambert models only."""
    if m.apo1 is not None and m.apo1 < m.alt1 or m.apo2 is not None and m.apo2 < m.alt2:
        raise ValueError("apoapsis altitude below periapsis altitude")
    di, dw, w1 = math.radians(m.di_deg), math.radians(m.dw_deg), math.radians(m.w1_deg)
    opts = []
    res, why = apse_transfer(m.r1, m.ra1, m.r2, m.ra2, di, dw, w1)
    o = Option("apse", "Apse-to-apse (optimal Δi split)", res["dv"] if res else math.nan, res["tof"] if res else math.nan,
               res["burns"] if res else 0, m.isp_chem, COMPLEXITY["apse"])
    if res:
        o.vec_plan, o.y0, o.legs = res["plan"], res["y0"], res["legs"]
        o.params = {k: res[k] for k in ("split_at_dep", "depart_deg", "rA", "rB")}
    else:
        o.applicable, o.note = False, why
    opts.append(o)
    res = lambert_optimal(m.r1, m.ra1, m.r2, m.ra2, di, dw, m.tof_max_h * 3600.0, w1=w1)
    o = Option("lambert", "Lambert (minimum Δv)", res["dv"] if res else math.nan, res["tof"] if res else math.nan,
               2, m.isp_chem, COMPLEXITY["lambert"])
    if res:
        o.vec_plan, o.y0, o.legs = res["plan"], res["y0"], res["legs"]
        o.params = {k: res[k] for k in ("th1_deg", "th2_deg", "time_limited")}
        o.note = "time limit sets the transfer time" if res["time_limited"] else ""
    else:
        o.applicable, o.note = False, "no single-revolution Lambert transfer within the time limit"
    opts.append(o)
    return opts


def _circular_options(m: Mission):
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

    dv, tof = low_thrust(r1, r2, di, m.accel_ep, m.isp_ep if m.mass0_kg is not None else None)
    o = Option("low_thrust", "Low-thrust spiral (Edelbaum)", dv, tof, 0, m.isp_ep, COMPLEXITY["low_thrust"])
    if di > EDELBAUM_MAX_DI:
        o.applicable, o.note = False, "Edelbaum's closed form holds only up to Δi = 114.6°"
    opts.append(o)
    return opts



# ------------------------------------------------- rule-based decision layer
def apply_rules(m: Mission, opts):
    tof_max = m.tof_max_h * 3600.0
    for o in opts:
        o.feasible, o.reasons = o.applicable, ([] if o.applicable else [o.note])
        if not o.applicable:
            continue
        need = (o.dv + o.dv_loss) * (1 + m.margin_pct / 100)
        if need > m.dv_budget:
            o.feasible = False
            o.reasons.append(f"Δv {o.dv:.3f} + loss {o.dv_loss:.3f} + {m.margin_pct:g}% margin = {need:.3f} "
                             f"> budget {m.dv_budget:.3f} km/s")
        if o.tof > tof_max * (1 + 1e-6):
            o.feasible = False
            o.reasons.append(f"TOF {o.tof/3600:.1f} h > limit {m.tof_max_h:.1f} h")
        if o.key == "low_thrust" and m.propulsion == "chemical":
            o.feasible = False
            o.reasons.append("no electric propulsion on board")
        if o.key != "low_thrust" and m.propulsion == "electric":
            o.feasible = False
            o.reasons.append("no high-thrust (chemical) propulsion on board")
        if m.mass0_kg is not None:          # the propellant on board must cover the same dv as the budget
            onboard = m.prop_ep_kg if o.key == "low_thrust" else m.prop_chem_kg
            prop = m.mass0_kg * (1 - math.exp(-need / (o.isp * G0)))
            if prop > onboard:
                o.feasible = False
                o.reasons.append(f"needs {prop:.0f} kg of propellant with the {m.margin_pct:g}% margin, "
                                 f"{onboard:.0f} kg on board")


def rule_advice(m: Mission):
    """Textbook rules of thumb used as the 'rule-based decision matrix'."""
    if m.elliptical:
        adv = ["Elliptical orbits: Hohmann, bi-elliptic, one-tangent and Edelbaum low thrust assume circular "
               "orbits, so the apse-to-apse and Lambert transfers are sized instead."]
        res, why = apse_transfer(m.r1, m.ra1, m.r2, m.ra2, *map(math.radians, (m.di_deg, m.dw_deg, m.w1_deg)))
        if res:
            adv.append("Apse lines aligned: a transfer between apsides is available, usually the cheapest.")
        else:
            adv.append(f"No apse-to-apse transfer: {why}. Only a Lambert transfer is offered.")
        return adv
    r1, r2 = m.r1, m.r2
    R = max(r1, r2) / min(r1, r2)
    adv = []
    if R < 11.94:
        adv.append(f"R = {R:.2f} < 11.94: without a plane change, Hohmann beats every bi-elliptic transfer.")
    elif R < 15.58:
        adv.append(f"11.94 ≤ R = {R:.2f} < 15.58: without a plane change, bi-elliptic wins only if rB is large enough.")
    else:
        adv.append(f"R = {R:.2f} ≥ 15.58: without a plane change, bi-elliptic beats Hohmann for any rB > r2.")
    if m.di_deg > 0:
        adv.append("Plane change: perform it where speed is lowest (largest radius) or combine with a burn. "
                   "That can make bi-elliptic the cheapest at any R.")
    if R < 1.0 + 1e-9 and m.di_deg > 48.94:
        adv.append("Δi > 48.94°: a three-impulse (bi-elliptic) plane change beats a single burn.")
    tof_h = hohmann(r1, r2, math.radians(m.di_deg))[1] / 3600
    if m.tof_max_h < tof_h:
        adv.append(f"TOF limit {m.tof_max_h:.2f} h < Hohmann {tof_h:.2f} h: a fast (one-tangent) transfer is required.")
    if m.propulsion != "chemical" and math.radians(m.di_deg) <= EDELBAUM_MAX_DI:
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


def _circular_burn(r, v, vafter, dinc, toward=None):
    """Velocity after a circular-model burn: speed vafter perpendicular to r, plane turned by
    dinc about r toward the plane with normal `toward`, or, without one, in the sense that
    increases inclination w.r.t. the initial plane (keeps split plane changes additive)."""
    h = np.cross(r, v)
    nh = h / np.linalg.norm(h)
    key = (lambda n: n[2]) if toward is None else (lambda n: -(n @ toward))
    n_new = min([_rotate(nh, r, dinc), _rotate(nh, r, -dinc)], key=key)
    return vafter * np.cross(n_new, r / np.linalg.norm(r))


def kepler(y, dt):
    """Two-body state dt seconds after state y: the universal-variable Kepler equation solved
    by Newton, then Lagrange f and g (Curtis, Algorithms 3.3 and 3.4)."""
    r0v, v0v = np.asarray(y[:3], float), np.asarray(y[3:], float)
    r0, sq = np.linalg.norm(r0v), math.sqrt(MU)
    vr0, alpha = np.dot(r0v, v0v) / r0, 2 / r0 - np.dot(v0v, v0v) / MU
    x = sq * abs(alpha) * dt
    for _ in range(100):
        z = alpha * x * x
        C, S = stumpff(z)
        F = r0 * vr0 / sq * x * x * C + (1 - alpha * r0) * x ** 3 * S + r0 * x - sq * dt
        dx = F / (r0 * vr0 / sq * x * (1 - z * S) + (1 - alpha * r0) * x * x * C + r0)
        x -= dx
        if abs(dx) < 1e-12 * max(1.0, abs(x)):
            break
    else:
        raise RuntimeError("kepler: Newton did not converge")
    z = alpha * x * x
    C, S = stumpff(z)
    rv = (1 - x * x / r0 * C) * r0v + (dt - x ** 3 * S / sq) * v0v
    r = np.linalg.norm(rv)
    return np.concatenate([rv, sq / (r * r0) * x * (z * S - 1) * r0v + (1 - x * x / r * C) * v0v])


def circular_vec_plan(r1, plan, di):
    """A circular model's burn plan as (y0, [(t, dv_vec)]) in the relative frame: the burns of
    simulate(), with the coasts solved analytically (kepler) instead of integrated, so the plan
    is exact. Plane changes turn toward the target plane, normal (0, -sin di, cos di), and must
    happen on the node line (x): if the first one wouldn't (a raising one-tangent transfer
    changes plane on arrival), the start point is moved back so that it does."""
    y, t0, phi = np.array([r1, 0, 0, 0, vcirc(r1), 0.0]), 0.0, 0.0
    for tb, vafter, dinc, *_ in plan:                  # in the start plane up to the first plane change
        if tb > t0:
            y, t0 = kepler(y, tb - t0), tb
        if dinc:
            phi = 0.0 if _on_node_line(math.atan2(y[1], y[0])) else math.atan2(y[1], y[0])
            break
        y[3:] = vafter * np.cross([0.0, 0.0, 1.0], y[:3] / np.linalg.norm(y[:3]))
    c, s, v1 = math.cos(phi), math.sin(phi), vcirc(r1)
    y, toward = np.array([r1 * c, -r1 * s, 0, v1 * s, v1 * c, 0.0]), np.array([0.0, -math.sin(di), math.cos(di)])
    y0, t0, out = y.copy(), 0.0, []
    for tb, vafter, dinc, *_ in plan:
        if tb > t0:
            y = kepler(y, tb - t0)
        vnew = _circular_burn(y[:3], y[3:], vafter, dinc, toward)
        out.append((tb, vnew - y[3:]))
        y, t0 = np.concatenate([y[:3], vnew]), tb
    return y0, out


def simulate(r1, plan, coast_after=0.0, n_out=2000):
    """Propagate an impulsive burn plan from a circular equatorial orbit at r1.
    Every burn in the ODST plans leaves velocity perpendicular to r (apsis or
    circularisation), so a burn = (time, speed after, plane rotation about r_hat)."""
    y = np.array([r1, 0, 0, 0, vcirc(r1), 0.0])
    t0, ts, ys, dv_used = 0.0, [], [], 0.0
    times = [b[0] for b in plan] + [plan[-1][0] + coast_after]
    for i, (tb, vafter, dinc, *_) in enumerate(plan):
        if tb > t0:
            sol = solve_ivp(_twobody, (t0, tb), y, rtol=1e-11, atol=1e-11, dense_output=True)
            tt = np.linspace(t0, tb, max(int(n_out * (tb - t0) / times[-1]), 50))
            ts.append(tt); ys.append(sol.sol(tt).T)
            y = sol.y[:, -1].copy()
        r, v = y[:3], y[3:]
        vnew = _circular_burn(r, v, vafter, dinc)
        dv_used += np.linalg.norm(vnew - v)
        y = np.concatenate([r, vnew])
        t0 = tb
    if coast_after > 0:
        sol = solve_ivp(_twobody, (t0, t0 + coast_after), y, rtol=1e-11, atol=1e-11, dense_output=True)
        tt = np.linspace(t0, t0 + coast_after, 400)
        ts.append(tt); ys.append(sol.sol(tt).T)
        y = sol.y[:, -1]
    return np.concatenate(ts), np.vstack(ys), y, dv_used


def simulate_vec(y0, plan, coast_after=0.0, n_out=2000):
    """Propagate a vector burn plan [(t, dv_vec)] from state y0 (two-body)."""
    y, t0, ts, ys, dv_used = np.array(y0, float), 0.0, [], [], 0.0
    t_end = (plan[-1][0] if plan else 0.0) + coast_after
    for tb, dv in plan + [(t_end, None)]:
        if tb > t0:
            sol = solve_ivp(_twobody, (t0, tb), y, rtol=1e-11, atol=1e-11, dense_output=True)
            tt = np.linspace(t0, tb, max(int(n_out * (tb - t0) / t_end), 50))
            ts.append(tt); ys.append(sol.sol(tt).T)
            y = sol.y[:, -1].copy()
        if dv is not None:
            y[3:] += dv
            dv_used += np.linalg.norm(dv)
        t0 = tb
    if not ts:
        ts, ys = [np.array([0.0])], [y[None, :]]
    return np.concatenate(ts), np.vstack(ys), y, dv_used


def target_mismatch(y, rp, ra, di, dw):
    """Relative error of angular momentum and absolute error of the eccentricity vector
    between state y and the orbit (rp, ra) in the plane turned by di, periapsis at dw."""
    r, v = y[:3], y[3:]
    h = np.cross(r, v)
    e = np.cross(v, h) / MU - r / np.linalg.norm(r)
    P, Q = orbit_basis(di, dw)
    ht = math.sqrt(MU * 2 * rp * ra / (rp + ra)) * np.cross(P, Q)
    return np.linalg.norm(h - ht) / np.linalg.norm(ht), np.linalg.norm(e - (ra - rp) / (ra + rp) * P)


def orbit_elements(y):
    r, v = y[:3], y[3:]
    h = np.cross(r, v)
    ev = np.cross(v, h) / MU - r / np.linalg.norm(r)
    a = 1.0 / (2 / np.linalg.norm(r) - np.dot(v, v) / MU)
    inc = math.degrees(math.acos(h[2] / np.linalg.norm(h)))
    return a, np.linalg.norm(ev), inc


# ------------------------------------------------------------------- export
GMAT_EPOCH = "01 Jan 2026 12:00:00.000"   # UTC; arbitrary, as the tool has no timing model yet


def _frame(m):
    return ("EarthMJ2000Eq, orbits from the elements entered" if m.R is not None else
            "EarthMJ2000Eq, the tool's frame: start orbit in the X-Y plane, node line on +X")


def inertial_plan(m, o):
    """Option o's start state and burns in EarthMJ2000Eq: the relative frame turned by m.R."""
    R = np.eye(3) if m.R is None else m.R
    return (np.concatenate([R @ o.y0[:3], R @ o.y0[3:]]), [(t, R @ dv) for t, dv in o.vec_plan])


def _num(x):
    return float(x) + 0.0               # plain float; + 0.0 turns -0.0 into 0.0


def _alts(m):
    return (m.alt1, m.alt1 if m.apo1 is None else m.apo1), (m.alt2, m.alt2 if m.apo2 is None else m.apo2)


def gmat_script(m: Mission, o: Option, report="odst_arrival.txt"):
    """NASA GMAT script that flies option o's impulsive burn plan from the tool's start state
    and reports the orbit it arrives on (same physics as the tool: point-mass Earth)."""
    g = lambda x: repr(_num(x))
    (p1, a1), (p2, a2) = _alts(m)
    y0, plan = inertial_plan(m, o)
    if m.R is None:
        orbits = [f"% Start orbit {p1:g} x {a1:g} km altitude; target {p2:g} x {a2:g} km, di {m.di_deg:g} deg, dw {m.dw_deg:g} deg"]
    else:
        e = m.elements
        orbits = [f"% Start orbit {p1:g} x {a1:g} km altitude, i {e[0]:g} deg, RAAN {e[1]:g} deg, argp {e[2]:g} deg",
                  f"% Target orbit {p2:g} x {a2:g} km altitude, i {e[3]:g} deg, RAAN {e[4]:g} deg, argp {e[5]:g} deg; "
                  f"planes {m.di_deg:.4g} deg apart"]
    s = ["% Orbital Manoeuvre Decision Support Tool: GMAT export",
         f"% Transfer: {o.name.replace('Δ', 'd')}, {len(plan)} burn(s), total dv {o.dv:.4f} km/s "
         f"(impulsive), flight time {o.tof / 3600:.2f} h", *orbits,
         f"% Frame: {_frame(m)}.",
         "%   The epoch is arbitrary: the tool has no timing model yet, so Sat starts at the departure point.",
         "% Physics as in the tool: point-mass Earth with the tool's mu and impulsive burns, so the",
         "%   Arrival report matches the target orbit to rounding. Edit TwoBody to add J2, drag or third bodies.",
         "", "% GMAT's default is 398600.4415, enough to shift a 10-day transfer's arrival by 0.2 km",
         f"Earth.Mu = {MU!r};", "", "Create Spacecraft Sat;", "Sat.DateFormat = UTCGregorian;", f"Sat.Epoch = '{GMAT_EPOCH}';",
         "Sat.CoordinateSystem = EarthMJ2000Eq;", "Sat.DisplayStateType = Cartesian;"]
    s += [f"Sat.{k} = {g(x)};" for k, x in zip(("X", "Y", "Z", "VX", "VY", "VZ"), y0)]
    sc = m.mass0_kg is not None
    if sc:                    # GMAT draws each burn's propellant from the tank: Sat.Fuel.FuelMass is reported
        s += ["", "% The spacecraft: the electric propellant is counted in the dry mass, as only chemical burns fly here",
              "Create ChemicalTank Fuel;", f"Fuel.FuelMass = {g(m.prop_chem_kg)};", "Fuel.AllowNegativeFuelMass = true;",
              "Fuel.FuelDensity = 1000;", f"Fuel.Volume = {g(m.prop_chem_kg / 1000)};",       # a tank just big enough
              f"Sat.DryMass = {g(m.mass0_kg - m.prop_chem_kg)};", "Sat.Tanks = {Fuel};"]
    s += ["", "Create ForceModel TwoBody;", "TwoBody.CentralBody = Earth;", "TwoBody.PointMasses = {Earth};",
          "Create Propagator Prop;", "Prop.FM = TwoBody;", "Prop.Type = RungeKutta89;", "Prop.InitialStepSize = 60;",
          "Prop.Accuracy = 1e-13;", "Prop.MinStep = 0;", "Prop.MaxStep = 600;"]
    for k, (_, dv) in enumerate(plan, 1):
        s += ["", f"Create ImpulsiveBurn Burn{k};", f"Burn{k}.CoordinateSystem = EarthMJ2000Eq;"]
        s += [f"Burn{k}.Element{j + 1} = {g(dv[j])};" for j in range(3)]
        if sc:
            s += [f"Burn{k}.DecrementMass = true;", f"Burn{k}.Tank = {{Fuel}};", f"Burn{k}.Isp = {g(m.isp_chem)};",
                  f"Burn{k}.GravitationalAccel = {g(G0 * 1e3)};"]
    s += ["", "Create ReportFile Arrival;", f"Arrival.Filename = '{report}';", "Arrival.Precision = 16;",
          "Arrival.WriteHeaders = true;", "", "Create OrbitView View;", "View.Add = {Sat, Earth};",
          "", "BeginMissionSequence;"]
    t = 0.0
    for k, (tb, _) in enumerate(plan, 1):
        if tb > t:
            s.append(f"Propagate Prop(Sat) {{Sat.ElapsedSecs = {g(tb - t)}}};")
        s.append(f"Maneuver Burn{k}(Sat);")
        t = tb
    s += ["Report Arrival " + " ".join(f"Sat.EarthMJ2000Eq.{k}" for k in ("X", "Y", "Z", "VX", "VY", "VZ"))
          + " Sat.Earth.SMA Sat.Earth.ECC Sat.EarthMJ2000Eq.INC Sat.EarthMJ2000Eq.RAAN Sat.EarthMJ2000Eq.AOP"
          + (" Sat.Fuel.FuelMass;" if sc else ";"),
          f"Propagate Prop(Sat) {{Sat.ElapsedSecs = {g(2 * math.pi * math.sqrt(((m.r2 + m.ra2) / 2) ** 3 / MU))}}};"]
    return "\n".join(s) + "\n"


def burns_json(m: Mission, o: Option):
    """Option o's burn plan as JSON: the start state, then each burn's time and inertial dv."""
    (p1, a1), (p2, a2) = _alts(m)
    y0, plan = inertial_plan(m, o)
    return json.dumps({
        "tool": "Orbital Manoeuvre Decision Support Tool", "transfer": o.name, "dv_kms": o.dv, "tof_s": o.tof,
        "start_orbit_alt_km": [p1, a1], "target_orbit_alt_km": [p2, a2], "di_deg": m.di_deg, "dw_deg": m.dw_deg,
        "w1_deg": m.w1_deg, "elements_deg": None if m.elements is None else [_num(x) for x in m.elements],
        "spacecraft": None if m.mass0_kg is None else {
            "mass0_kg": m.mass0_kg, "thrust_chem_N": m.thrust_chem_N, "isp_chem_s": m.isp_chem, "prop_chem_kg": m.prop_chem_kg,
            "propellant_kg": _num(o.prop_kg)},
        "frame": _frame(m), "epoch_utc": GMAT_EPOCH, "mu_km3_s2": MU,
        "start_state_km_kms": [_num(x) for x in y0],
        "burns": [{"t_s": _num(t), "dv_kms": [_num(x) for x in dv], "dv_mag_kms": _num(np.linalg.norm(dv))}
                  for t, dv in plan]}, indent=2, ensure_ascii=False) + "\n"


# ------------------------------------------------------------------- report
def print_table(m: Mission, opts):
    o1 = f"{m.alt1:.0f}" if m.apo1 is None else f"{m.alt1:.0f} x {m.apo1:.0f}"
    o2 = f"{m.alt2:.0f}" if m.apo2 is None else f"{m.alt2:.0f} x {m.apo2:.0f}"
    dw = f", Δω = {m.dw_deg:g}°" if m.elliptical else ""
    print(f"\nMission: alt {o1} -> {o2} km, Δi = {m.di_deg:.2f}°{dw}, "
          f"budget {m.dv_budget} km/s (margin {m.margin_pct:g}%, T/W {m.tw_chem:g}), TOF max {m.tof_max_h:.1f} h, weights {m.weights}")
    hdr = f"{'Option':34s}{'Δv km/s':>9s}{'loss m/s':>9s}{'TOF':>12s}{'burns':>6s}{'mp/m0':>8s}{'cplx':>6s}{'WS':>7s}{'TOPSIS':>8s}  status"
    print(hdr); print("-" * len(hdr))
    for o in opts:
        tof = f"{o.tof/3600:.2f} h" if o.tof < 72 * 3600 else f"{o.tof/86400:.1f} d"
        st = "OK" if o.feasible else "X " + "; ".join(o.reasons)
        ws = f"{o.score_ws:.3f}" if o.feasible else "-"
        tp = f"{o.score_topsis:.3f}" if o.feasible else "-"
        print(f"{o.name:34s}{o.dv:9.4f}{1e3 * o.dv_loss:9.1f}{tof:>12s}{o.burns:6d}{o.prop_frac:8.3f}{o.complexity:6.1f}{ws:>7s}{tp:>8s}  {st}")
    best = recommend(opts)
    print("Recommended:", best.name if best else "NONE FEASIBLE - relax constraints")
    for a in rule_advice(m):
        print("  rule:", a)


if __name__ == "__main__":
    print_table(Mission(), evaluate(Mission()))
