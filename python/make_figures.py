"""Generate the report figures into ../docs/figures.

Uses the shared ODST design tokens (same as web/index.html, matlab/odst_theme.m
and the report): white canvas, left-aligned near-black titles, grey axis text,
hairline horizontal grid, no box or tick marks, one fixed colour per manoeuvre.
"""
import math
import os
import sys
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyBboxPatch, Patch

from odst import (RE, WEIGHT_PRESETS, Mission, bielliptic, evaluate, hohmann,
                  recommend, simulate, vcirc)

sys.stdout.reconfigure(encoding="utf-8")
OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "figures")
os.makedirs(OUT, exist_ok=True)

# ---- design tokens (light) -------------------------------------------------
LABEL, LABEL2, LABEL3 = "#1d1d1f", "#6e6e73", "#86868b"
GRID, AXIS, SURF, PAGE = "#ececf0", "#d2d2d7", "#ffffff", "#f5f5f7"
EARTH, ORBIT, NONE = "#a9cbf3", "#aeaeb2", "#e5e5ea"
ACCENT, ACCENT_SOFT, GOOD_SOFT, BAD_SOFT, WARN_SOFT = "#0071e3", "#e8f2fd", "#e7f2ea", "#fcebec", "#fff0db"
# Entity -> colour (fixed, never by rank)
COL = {"hohmann": "#2a78d6", "hohmann_split": "#eb6834", "bielliptic": "#1baf7a",
       "one_tangent": "#eda100", "low_thrust": "#e87ba4"}
BI_RAMP = ["#46c094", "#1baf7a", "#13875e", "#0b5c40"]   # ordinal ramp of the bi-elliptic hue
SHORT = {"hohmann": "Hohmann", "hohmann_split": "Hohmann\nsplit", "bielliptic": "Bi-\nelliptic",
         "one_tangent": "One-\ntangent", "low_thrust": "Low-\nthrust"}


def label(k):
    return SHORT[k].replace("-\n", "-").replace("\n", " ")


plt.rcParams.update({
    "font.family": ["Segoe UI", "Arial", "DejaVu Sans"], "font.size": 10,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.labelcolor": LABEL2, "axes.labelsize": 10,
    "axes.titlecolor": LABEL, "axes.titlesize": 12, "axes.titleweight": "semibold",
    "axes.titlelocation": "left", "axes.titlepad": 10,
    "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
    "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.8,
    "xtick.color": LABEL2, "ytick.color": LABEL2, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "xtick.major.size": 0, "ytick.major.size": 0, "xtick.minor.size": 0, "ytick.minor.size": 0,
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "legend.frameon": False, "legend.fontsize": 9, "figure.dpi": 150,
})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print("wrote", name)


def heading(fig, title, sub=None, y=1.0):
    """Left-aligned figure title + optional grey subtitle."""
    fig.text(0.01, y, title, ha="left", va="bottom", fontsize=14, fontweight="bold", color=LABEL)
    if sub:
        fig.text(0.01, y - 0.012, sub, ha="left", va="top", fontsize=10.5, color=LABEL2)


def title_sub(ax, title, sub):
    """Axes title with a grey subtitle line between it and the plot."""
    ax.set_title(title, pad=26)
    ax.annotate(sub, xy=(0, 1), xycoords="axes fraction", xytext=(0, 8), textcoords="offset points",
                va="bottom", ha="left", fontsize=9.5, color=LABEL2)


def fmt_t(s):
    return f"{s/3600:.2f} h" if s < 72 * 3600 else f"{s/86400:.1f} days"


# ---------------------------------------------------------------- Fig 1
def fig_case_bars(m, name, title, sub):
    opts = evaluate(m)
    rec = recommend(opts)
    fig, axs = plt.subplots(1, 4, figsize=(13.6, 3.4))
    panels = [("Total Δv", "km/s", lambda o: o.dv, False),
              ("Propellant", "fraction of initial mass", lambda o: o.prop_frac, False),
              ("Transfer time", "hours, log scale", lambda o: o.tof / 3600, True),
              ("Score", "weighted sum, 1 = best", lambda o: o.score_ws if o.feasible else 0, False)]
    labels = [SHORT[o.key] for o in opts]
    for ax, (t, unit, f, logy) in zip(axs, panels):
        vals = [f(o) for o in opts]
        bars = ax.bar(range(len(opts)), vals, width=0.62, color=[COL[o.key] for o in opts],
                      edgecolor=SURF, linewidth=2, zorder=3)
        for b, o in zip(bars, opts):
            if not o.feasible:
                b.set_alpha(0.25)
        ax.set_title(t)
        ax.set_ylabel(unit, fontsize=9, color=LABEL3)
        ax.set_xticks(range(len(opts)), labels, fontsize=8.2)
        ax.spines["bottom"].set_color(AXIS)
        if logy:
            ax.set_yscale("log")
            ax.minorticks_off()
            ax.set_yticks([1, 10, 100, 1000], ["1", "10", "100", "1000"])
        if t == "Score":
            ax.set_ylim(0, 1.22)
            ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
            for i, o in enumerate(opts):
                if o is rec:
                    ax.annotate("Recommended", (i, o.score_ws), xytext=(0, 5), textcoords="offset points",
                                ha="center", fontsize=8.5, fontweight="semibold", color=LABEL)
    heading(fig, title, sub, y=1.07)
    fig.text(0.01, -0.05, "Faded bars are infeasible or not applicable under the mission constraints.",
             fontsize=9, color=LABEL3)
    fig.tight_layout()
    save(fig, name)
    return opts


# ---------------------------------------------------------------- Fig 2
def burn_points(r1, plan):
    pts = []
    for i in range(len(plan)):
        if plan[i][0] > 0:
            _, _, y, _ = simulate(r1, plan[: i + 1], coast_after=0.0)
        else:
            y = np.array([r1, 0, 0])
        pts.append(y[:2])
    return np.array(pts)


def fig_trajectories(m, name, title, sub, keys):
    opts = {o.key: o for o in evaluate(m)}
    fig, axs = plt.subplots(1, len(keys), figsize=(4.1 * len(keys), 4.2))
    th = np.linspace(0, 2 * np.pi, 400)
    k = 1e-3                                           # 10^3 km
    for ax, key in zip(np.atleast_1d(axs), keys):
        o = opts[key]
        _, Y, _, _ = simulate(m.r1, o.burn_plan, coast_after=0.0, n_out=3000)
        ax.fill(k * RE * np.cos(th), k * RE * np.sin(th), color=EARTH, zorder=2, lw=0)
        for r, ls in [(m.r1, "-"), (m.r2, "--")]:
            ax.plot(k * r * np.cos(th), k * r * np.sin(th), color=ORBIT, lw=1, ls=ls, zorder=2)
        ax.plot(k * Y[:, 0], k * Y[:, 1], color=COL[key], lw=2.2, zorder=4)
        pts = burn_points(m.r1, o.burn_plan) * k
        ax.scatter(pts[:, 0], pts[:, 1], s=46, color=COL[key], edgecolor=SURF, linewidth=1.8, zorder=5)
        ax.set_aspect("equal")
        ax.grid(True, axis="both")
        title_sub(ax, o.name, f"Δv = {o.dv:.3f} km/s · TOF = {fmt_t(o.tof)}")
        ax.set_xlabel("x [10³ km]"); ax.set_ylabel("y [10³ km]")
    heading(fig, title, sub, y=1.06)
    fig.tight_layout()
    save(fig, name)


# ---------------------------------------------------------------- Fig 3
def fig_hohmann_vs_bielliptic():
    r1 = 7000.0
    R = np.geomspace(1.5, 80, 300)
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    dh = np.array([hohmann(r1, x * r1, 0)[0] for x in R]) / vcirc(r1)
    ax.plot(R, dh, color=COL["hohmann"], lw=2.4, label="Hohmann", zorder=4)
    for c, n in zip(BI_RAMP, [2, 5, 15, 1e5]):
        db = np.array([bielliptic(r1, x * r1, n * x * r1, 0)[0] for x in R]) / vcirc(r1)
        lab = "Bi-elliptic, rB = ∞" if n > 1e4 else f"Bi-elliptic, rB = {n:g}·r2"
        ax.plot(R, db, color=c, lw=1.8, ls="--" if n > 1e4 else "-", label=lab)
    for x in (11.94, 15.58):
        ax.axvline(x, color=LABEL3, lw=0.9, ls=":")
        ax.text(x * 1.03, 0.405, f"R = {x}", color=LABEL2, fontsize=9, rotation=90, va="bottom")
    ax.set_xscale("log")
    ax.minorticks_off()
    ax.set_xticks([2, 5, 10, 20, 50], ["2", "5", "10", "20", "50"])
    ax.set_xlabel("Radius ratio R = r2 / r1"); ax.set_ylabel("Total Δv / circular speed at r1")
    ax.set_ylim(0.4, 0.58)
    title_sub(ax, "Hohmann vs bi-elliptic",
              "Coplanar transfers. Bi-elliptic wins above R = 15.58 for any rB, above 11.94 only for large rB.")
    ax.legend(loc="lower right")
    save(fig, "fig3_hohmann_vs_bielliptic.png")


# ---------------------------------------------------------------- Fig 4
def fig_plane_split(m):
    r1, r2, di = m.r1, m.r2, math.radians(m.di_deg)
    s = np.linspace(0, 1, 401)
    dv = np.array([hohmann(r1, r2, di, x)[0] for x in s])
    i = int(np.argmin(dv))
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    ax.plot(s * m.di_deg, dv, color=COL["hohmann_split"], lw=2.4)
    ax.scatter([s[i] * m.di_deg], [dv[i]], s=70, color=COL["hohmann_split"], edgecolor=SURF, lw=2, zorder=5)
    ax.annotate(f"Optimum: {s[i]*m.di_deg:.2f}° at the LEO burn\nΔv = {dv[i]:.4f} km/s, saves {1000*(dv[0]-dv[i]):.0f} m/s",
                (s[i] * m.di_deg, dv[i]), xytext=(28, 46), textcoords="offset points", fontsize=9.5, color=LABEL,
                arrowprops=dict(arrowstyle="-", color=LABEL3, lw=0.8))
    ax.set_xlabel("Plane change performed at the LEO (perigee) burn [°]")
    ax.set_ylabel("Total Δv [km/s]")
    title_sub(ax, f"Optimal plane-change split (Δi = {m.di_deg}°)",
              "Case A, Hohmann transfer LEO → GEO, minimised numerically")
    save(fig, "fig4_plane_change_split.png")


# ---------------------------------------------------------------- Fig 5
KEYS = ["hohmann", "hohmann_split", "bielliptic", "one_tangent", "low_thrust"]


def _rec_grid(xs, ys, make):
    Z = np.zeros((len(ys), len(xs)))
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            rec = recommend(evaluate(make(x, y)))
            Z[i, j] = KEYS.index(rec.key) if rec else len(KEYS)
    return Z


def fig_decision_map():
    alt1 = 400.0
    r1 = RE + alt1
    Rs = np.geomspace(1.05, 50, 44)
    cmap = ListedColormap([COL[k] for k in KEYS] + [NONE])
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    tofs = np.geomspace(0.3, 24 * 365, 40)
    Za = _rec_grid(Rs, tofs, lambda R, t: Mission(alt1=alt1, alt2=R * r1 - RE, di_deg=28.5, dv_budget=10,
                   tof_max_h=t, propulsion="both", weights=WEIGHT_PRESETS["Propellant-critical"]))
    dis = np.linspace(0, 90, 37)
    Zb = _rec_grid(Rs, dis, lambda R, d: Mission(alt1=alt1, alt2=R * r1 - RE, di_deg=d, dv_budget=12,
                   tof_max_h=24 * 365, propulsion="chemical", weights=WEIGHT_PRESETS["Propellant-critical"]))
    for ax, Z, ys, yl, t in [
            (axs[0], Za, tofs, "Transfer-time limit [h]", "(a) Δi = 28.5°, chemical + electric"),
            (axs[1], Zb, dis, "Inclination change Δi [°]", "(b) Time limit 1 year, chemical only")]:
        ax.pcolormesh(Rs, ys, Z, cmap=cmap, vmin=-0.5, vmax=len(KEYS) + 0.5, shading="nearest",
                      edgecolors="none", linewidth=0)
        ax.set_xscale("log")
        ax.minorticks_off()
        ax.set_xticks([1.05, 2, 5, 10, 20, 50], ["1", "2", "5", "10", "20", "50"])
        ax.set_xlabel("Radius ratio R = r2 / r1  (r1 = LEO 400 km)")
        ax.set_ylabel(yl)
        ax.set_title(t)
        ax.grid(False)
    axs[0].set_yscale("log")
    axs[0].minorticks_off()
    axs[0].set_yticks([1, 10, 100, 1000], ["1 h", "10 h", "100 h", "1000 h"])
    axs[0].axvline(6.22, color=LABEL, lw=0.9, ls=":")
    axs[0].text(6.22 * 1.04, 0.36, "GEO", fontsize=9, color=LABEL, fontweight="semibold")
    handles = [Patch(color=COL[k], label=label(k)) for k in KEYS] + [Patch(color=NONE, label="None feasible")]
    fig.legend(handles=handles, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.06), fontsize=9.5)
    heading(fig, "Decision maps: recommended manoeuvre across the mission space",
            "Propellant-critical priorities. Each cell is one full evaluation of all five manoeuvres.", y=1.05)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, "fig5_decision_map.png")

    # single-panel versions for the report layout (panel a pairs with MATLAB's map)
    for tag, Z, ys, yl, t, logy in [
            ("a", Za, tofs, "Transfer-time limit [h]", "Δi = 28.5°, chemical + electric", True),
            ("b", Zb, dis, "Inclination change Δi [°]", "Time limit 1 year, chemical only", False)]:
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
        ax.pcolormesh(Rs, ys, Z, cmap=cmap, vmin=-0.5, vmax=len(KEYS) + 0.5, shading="nearest",
                      edgecolors="none", linewidth=0)
        ax.set_xscale("log")
        ax.minorticks_off()
        ax.set_xticks([1.05, 2, 5, 10, 20, 50], ["1", "2", "5", "10", "20", "50"])
        ax.set_xlabel("Radius ratio R = r2 / r1  (r1 = LEO 400 km)")
        ax.set_ylabel(yl)
        ax.grid(False)
        if logy:
            ax.set_yscale("log")
            ax.minorticks_off()
            ax.set_yticks([1, 10, 100, 1000], ["1 h", "10 h", "100 h", "1000 h"])
        title_sub(ax, "Recommended manoeuvre", t + ", propellant-critical")
        present = sorted(set(int(v) for v in np.unique(Z)))
        hs = [Patch(color=COL[KEYS[i]], label=label(KEYS[i])) if i < len(KEYS) else Patch(color=NONE, label="None feasible")
              for i in present]
        ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=9)
        save(fig, f"fig5{tag}_decision_map.png")


# ---------------------------------------------------------------- Fig 6
def fig_weight_sensitivity(m):
    base = np.array(WEIGHT_PRESETS["Balanced"], float)
    others = np.delete(base, 2) / np.delete(base, 2).sum()
    wt = np.linspace(0, 0.9, 91)
    curves = {}
    for w in wt:
        weights = np.insert(others * (1 - w), 2, w)
        for o in evaluate(replace(m, weights=tuple(weights))):
            if o.feasible:
                curves.setdefault(o.key, []).append((w, o.score_ws))
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    for k, pts in curves.items():
        pts = np.array(pts)
        ax.plot(pts[:, 0], pts[:, 1], color=COL[k], lw=2.2, label=label(k))
    ax.axvline(0.863, color=LABEL3, lw=0.9, ls=":")
    ax.text(0.855, 0.05, "0.86", color=LABEL2, fontsize=9, ha="right")
    ax.set_xlim(0, 0.92)
    ax.set_xlabel("Weight on transfer time (other weights scaled proportionally)")
    ax.set_ylabel("Weighted-sum score")
    title_sub(ax, "Sensitivity to the transfer-time weight",
              "Case A, LEO → GEO with Δi = 28.5°: one-tangent overtakes Hohmann above 0.86")
    ax.legend(loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.32))
    save(fig, "fig6_weight_sensitivity.png")


# ---------------------------------------------------------------- Fig 7 flowchart
def fig_flowchart():
    fig, ax = plt.subplots(figsize=(10, 7.6))
    ax.set_xlim(0, 100); ax.set_ylim(0, 107); ax.axis("off")

    def box(x, y, w, h, head, body, fc=PAGE):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.35,rounding_size=1.6",
                                    fc=fc, ec="none"))
        ax.text(x, y + h / 2 - 2.0 if body else y, head, ha="center", va="center", fontsize=9.5,
                fontweight="semibold", color=LABEL)
        if body:
            ax.text(x, y + h / 2 - 3.8, body, ha="center", va="top", fontsize=8.6, color=LABEL2, linespacing=1.35)

    def arrow(x1, y1, x2, y2):
        ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(arrowstyle="-|>", color=LABEL3, lw=1))

    box(50, 94, 70, 7, "Mission inputs", "r1, r2, Δi · Δv budget · time limit · propulsion · Isp, F/m · priority weights")
    box(50, 82, 78, 7, "Model five candidates", "Hohmann · Hohmann with Δi split · bi-elliptic · one-tangent · Edelbaum low thrust")
    box(50, 70, 78, 7, "Optimise each one", "split s* = argmin Δv · rB* = argmin Δv with TOF ≤ limit · k* = smallest k meeting the limit")
    box(50, 58, 70, 7, "Rule filters", "Δv ≤ budget · TOF ≤ limit · propulsion on board · option applicable")
    box(21, 44, 36, 11, "Advisory rules",
        "R < 11.94: Hohmann optimal\nR > 15.58: bi-elliptic cheaper\nplane change at the largest radius", fc=ACCENT_SOFT)
    box(79, 44, 36, 11, "Time rule",
        "time limit below Hohmann TOF:\nfast one-tangent transfer\nlow thrust if time allows", fc=ACCENT_SOFT)
    box(50, 29.5, 70, 9, "Score the feasible options",
        "weighted sum (min–max) with TOPSIS cross-check\non Δv, propellant, log TOF, complexity")
    box(50, 17, 40, 5.6, "Any feasible option?", "", fc=WARN_SOFT)
    box(22, 5, 38, 7, "No", "report the violated constraints", fc=BAD_SOFT)
    box(78, 5, 38, 7, "Yes", "recommend the top score, show ranking and plots", fc=GOOD_SOFT)
    arrow(50, 90.3, 50, 85.7); arrow(50, 78.3, 50, 73.7); arrow(50, 66.3, 50, 61.7)
    arrow(40, 54.3, 26, 49.7); arrow(60, 54.3, 74, 49.7); arrow(50, 54.3, 50, 34.2)
    arrow(26, 38.3, 40, 34.2); arrow(74, 38.3, 60, 34.2); arrow(50, 24.8, 50, 20)
    arrow(40, 14, 28, 8.7); arrow(60, 14, 72, 8.7)
    ax.text(0, 105.5, "Decision logic", fontsize=14, fontweight="bold", color=LABEL, va="bottom")
    ax.text(0, 104.6, "Rule-based matrix, numerical optimisation and multi-criteria evaluation in one pipeline",
            fontsize=10.5, color=LABEL2, va="top")
    save(fig, "fig7_flowchart.png")


if __name__ == "__main__":
    leo_geo = Mission(alt1=400, alt2=35786, di_deg=28.5, dv_budget=6.5, tof_max_h=24 * 200)
    high_R = Mission(alt1=7000 - RE, alt2=105000 - RE, di_deg=0, dv_budget=6.0, tof_max_h=24 * 60,
                     propulsion="chemical", weights=WEIGHT_PRESETS["Propellant-critical"])
    rapid = Mission(alt1=400, alt2=800, di_deg=0, dv_budget=1.0, tof_max_h=0.6,
                    weights=WEIGHT_PRESETS["Time-critical"])
    fig_case_bars(leo_geo, "fig1a_case_leo_geo.png", "Case A: LEO 400 km → GEO, Δi = 28.5°",
                  "Balanced priorities, chemical and electric propulsion available")
    fig_case_bars(high_R, "fig1b_case_high_ratio.png", "Case B: 7 000 km → 105 000 km (R = 15)",
                  "Propellant-critical priorities, chemical propulsion only")
    fig_case_bars(rapid, "fig1c_case_rapid.png", "Case C: LEO 400 km → 800 km in 36 minutes",
                  "Time-critical priorities, Δv budget 1 km/s")
    fig_trajectories(high_R, "fig2_trajectories_case_b.png", "Case B transfers, numerically propagated",
                     "Two-body motion with impulsive burns; markers show each burn", ["hohmann", "bielliptic", "one_tangent"])
    fig_hohmann_vs_bielliptic()
    fig_plane_split(leo_geo)
    fig_decision_map()
    fig_weight_sensitivity(leo_geo)
    fig_flowchart()
