# Orbital Manoeuvre Decision Support Tool

[![tests](https://github.com/satwiksharma01/orbital-manoeuvre-decision-support/actions/workflows/tests.yml/badge.svg)](https://github.com/satwiksharma01/orbital-manoeuvre-decision-support/actions/workflows/tests.yml)

**Which manoeuvre should a mission fly?** A decision-support tool that sizes five orbital transfers, optimises each one numerically, removes those that break the mission's constraints and ranks the rest by the mission's priorities. Built in MATLAB and Simulink, cross-checked in Python, with an interactive web planner.

![Web planner demo: picking example missions re-ranks the manoeuvres and draws each transfer](assets/demo.gif)

**[Try the live demo](https://satwiksharma01.github.io/orbital-manoeuvre-decision-support/web/)**

## What it does

Given an initial and a target circular orbit, an inclination change, a Δv budget, a time limit, the propulsion on board and the mission's priorities, the tool:

1. **Models five candidates:** Hohmann, Hohmann with an optimised plane-change split, bi-elliptic, one-tangent burn and an Edelbaum low-thrust spiral.
2. **Optimises each one:** the plane-change split (`fminbnd`), the bi-elliptic apoapsis under a time limit (grid search + `fminbnd` + `fzero`) and the cheapest fast transfer that meets a deadline (`fzero`).
3. **Applies rules:** hard constraints (Δv with a margin, default 5 %, time, propulsion) plus textbook rules such as the 11.94 / 15.58 Hohmann vs bi-elliptic thresholds.
4. **Ranks the survivors:** weighted-sum scoring cross-checked with TOPSIS, with presets for balanced, propellant-critical, time-critical and low-complexity missions.

## Results

| Mission | Recommendation | Why |
|---|---|---|
| LEO 400 km → GEO, Δi = 28.5° | Hohmann, 4.22 km/s, 5.3 h | Optimal split saves only 25 m/s; stays first until the time weight exceeds 0.86 |
| Same, propellant-critical, 1-year limit | Low-thrust spiral | Propellant fraction 0.31 instead of 0.74 |
| 7 000 → 105 000 km (R = 15) | Hohmann | Bi-elliptic saves 1.2 % Δv but takes 12 days instead of 18 hours |
| 400 → 800 km within 36 min | One-tangent, k* = 1.0098 | Hohmann takes 48 min; root-finding gives the cheapest transfer that arrives in time |

## Verification

- Textbook benchmarks (Hohmann LEO→GEO 3.893 km/s, the 11.94 bi-elliptic threshold, plane-change limits).
- Every burn plan propagated with `ode45`: final orbit within 0.5 km, e < 10⁻⁹, Δv matches the model to 10⁻⁵ km/s.
- **Simulink** finite-burn model: the error grows as τ², about 56 km at τ = 60 s for LEO→GEO, so the impulsive model holds for short chemical burns.
- An independent Python implementation and a Node test of the web engine agree with MATLAB to every printed digit.

| MATLAB GUI | Simulink model |
|---|---|
| ![MATLAB app](assets/matlab_app.png) | ![Simulink model](assets/simulink_model.png) |

![Case A results](assets/case_A_bars.png)
![Decision map](assets/decision_map.png)
![Propagated transfers](assets/case_B_trajectories.png)
![Simulink finite-burn verification](assets/simulink_verification.png)

## Scope and limitations

This is an early trade-study and teaching tool: it tells you which family of transfer to fly, roughly what Δv it costs and roughly how long it takes. It is not a mission-analysis or operations tool. It assumes:

| Assumption | What it leaves out |
|---|---|
| Start and target orbits are circular | Elliptical targets such as GTO or Molniya |
| Two-body gravity | J2 (RAAN drifts by degrees a day in LEO) and drag below about 500 km |
| Burns are instantaneous | Gravity losses. The Simulink check bounds the error, but none is added to the Δv budget |
| Burn timing is ignored | Phasing, rendezvous and launch windows |
| Low thrust is Edelbaum at constant acceleration | Power limits, eclipses and acceleration rising as propellant is used |
| A fixed Δv margin (default 5 %) | Uncertainty in Isp, thrust and pointing |
| Weights are chosen by the user | Only the time weight is swept; other weights can flip the ranking |

Every check so far compares the tool with textbooks and with itself, not yet with an industry tool such as GMAT.

## Run it

**MATLAB** (R2021a+, Simulink optional), from the `matlab` folder:

```matlab
odst_tests   % validation suite
odst_main    % case studies, figures (matlab/results) and the Simulink check
odst_app     % interactive planner
```

**Web:** open `web/index.html` in any browser. No build step and no dependencies.

**Python** (3.10+, `numpy scipy matplotlib`):

```bash
python python/validate.py        # cross-check
node web/test_engine.js          # web engine regression test
```

## Repository layout

| Path | Contents |
|---|---|
| `matlab/` | Models, optimisers, decision pipeline (`odst_evaluate.m`), GUI, Simulink model and builder, tests |
| `python/` | Reference implementation, validation and figure scripts |
| `web/` | Self-contained web planner and its engine test |
| `assets/` | Screenshots and result figures used in this README |
