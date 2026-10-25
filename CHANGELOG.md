# Changelog

All notable changes to this project. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed

- **Result figures:** the case A and B bars were regenerated, so their propellant fractions include the gravity loss (finding 9). The decision map came out identical.
- **Repository:** removed an unused screenshot and the backup copy Simulink kept when it upgraded the model; `*.slx.r20*` is now ignored.

## [1.0.0] - 2026-10-24

The first tagged release: the decision-support tool in MATLAB/Simulink, Python and the web, for circular and elliptical orbits, checked against textbooks, numerical propagation, NASA GMAT and each other.

### Added

- **Five transfers between circular orbits:** Hohmann, Hohmann with an optimised plane-change split, bi-elliptic with an optimised apoapsis, one-tangent burn (the cheapest that meets a deadline) and an Edelbaum low-thrust spiral.
- **Elliptical orbits:** periapsis and apoapsis for each orbit and the target's apse angle Δω. Two transfers for elliptical missions: apse-to-apse with an optimised plane-change split, and the minimum-Δv single-revolution Lambert transfer (universal-variable solver, grid search, golden section, Nelder–Mead).
- **Decision layer:** hard rules (Δv budget, time limit, propulsion on board), textbook rules, weighted-sum ranking cross-checked with TOPSIS, and four priority presets.
- **Δv margin** (default 5 %) and a **finite-burn gravity-loss estimate** (Robbins) from a thrust-to-weight input; the budget check covers (Δv + loss) × (1 + margin).
- **Lambert ground clearance:** every Lambert arc stays at least 100 km above the surface, enforced as an exact penalty.
- **Three implementations that agree:** MATLAB (engine, app, Simulink finite-burn model), Python (reference engine, validation, figures) and a self-contained web planner, live on GitHub Pages.
- **Verification:**
  - `python/validate.py` (52 checks), `matlab/odst_tests.m` (64) and `web/test_engine.js` (59), all in CI on every push.
  - A NASA GMAT cross-check (`python/gmat_crosscheck.py`, `gmat/`), run locally.
- **README:** results, verification, scope and limitations, an animated demo of the web planner and `tools/record_demo_gif.py` to re-record it.

### Fixed

Each fix below is in the release; the numbers are the project's internal finding numbers.

- **Split plane change:** its second half rotated the wrong way in the propagated burn plan (finding 1).
- **Exact deadline:** a transfer solved to meet the time limit exactly was rejected for missing it (finding 2).
- **Web form:** valid inputs were marked invalid (finding 3).
- **Dark figures:** MATLAB R2025a and later drew every figure dark (finding 4).
- **Phone layout:** the web planner had no doctype or viewport tag, so it rendered in quirks mode and phones got the desktop layout (finding 7).
- **Edelbaum limit:** the formula was applied beyond its 114.6° validity; low thrust is no longer offered above it (finding 8).
- **Propellant fraction:** it now includes the gravity loss, as the budget check already did (finding 9).
- **Rule text:** the 11.94 / 15.58 Hohmann vs bi-elliptic rules now say they hold for coplanar transfers (finding 10).
- **Ground clearance:** Lambert optima could pass through the Earth with low perigees and tight time limits (finding 11).

### Known limitations

See *Scope and limitations* in the README:
- two-body gravity;
- impulsive burns with an estimated loss;
- no phasing or waiting time;
- the plane change pivots on the start orbit's apse line;
- single-revolution Lambert, a grid-plus-refinement search (up to 0.35 % above the best transfer found when a time limit binds between low-perigee orbits);
- Edelbaum at constant acceleration, between circular orbits only.

[Unreleased]: https://github.com/satwiksharma01/orbital-manoeuvre-decision-support/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/satwiksharma01/orbital-manoeuvre-decision-support/releases/tag/v1.0.0
