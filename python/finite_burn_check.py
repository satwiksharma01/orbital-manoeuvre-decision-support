"""Replica of matlab/odst_sim_dynamics.m (finite-burn Simulink model) for verification."""
import math, sys
import numpy as np
from scipy.integrate import solve_ivp
from odst import MU, RE, Mission, evaluate, recommend, orbit_elements, _rotate
sys.stdout.reconfigure(encoding="utf-8")

def coast(r, v, dt, n=50):
    """Fixed-step RK4 two-body coast (mirrors the MATLAB helper)."""
    y = np.concatenate([r, v]); h = dt / n
    g = lambda y: np.concatenate([y[3:], -MU * y[:3] / np.linalg.norm(y[:3]) ** 3])
    for _ in range(n):
        k1 = g(y); k2 = g(y + h / 2 * k1); k3 = g(y + h / 2 * k2); k4 = g(y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return y[:3], y[3:]


def run(o, r1, r2, tau):
    plan = o.burn_plan
    acc = [None] * len(plan)
    def f(t, y):
        r, v = y[:3], y[3:]
        a = -MU * r / np.linalg.norm(r) ** 3
        for k, (tb, va, dinc, *_) in enumerate(plan):
            if tb - tau / 2 <= t < tb + tau / 2:
                if acc[k] is None:
                    rc, vc = coast(r, v, tb - t)   # state at the burn centre (no thrust)
                    nh = np.cross(rc, vc); nh /= np.linalg.norm(nh)
                    n = min([_rotate(nh, rc, dinc), _rotate(nh, rc, -dinc)], key=lambda n: n[2])
                    acc[k] = (va * np.cross(n, rc / np.linalg.norm(rc)) - vc) / tau
                a = a + acc[k]
        return np.concatenate([v, a])
    n1 = math.sqrt(MU / r1 ** 3); th = -n1 * tau; v1 = math.sqrt(MU / r1)
    x0 = [r1 * math.cos(th), r1 * math.sin(th), 0, -v1 * math.sin(th), v1 * math.cos(th), 0]
    tend = plan[-1][0] + tau + 2 * math.pi * math.sqrt(r2 ** 3 / MU)
    sol = solve_ivp(f, (-tau, tend), x0, rtol=1e-9, atol=1e-9, max_step=tau / 4)
    return orbit_elements(sol.y[:, -1])

if __name__ == "__main__":
    for m in [Mission(), Mission(alt1=7000 - RE, alt2=105000 - RE, di_deg=0)]:
        for o in evaluate(m):
            if o.burn_plan and o.applicable:
                for tau in (10, 60, 300):
                    a, e, i = run(o, m.r1, m.r2, tau)
                    print(f"{o.name:30s} tau={tau:4d}s  a err {a - m.r2:+10.2f} km  e {e:.1e}  i err {i - m.di_deg:+.3f} deg")
