# Changelog

All notable changes to this project. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **Spacecraft model:** a wet mass, chemical thrust, electric power and efficiency, and the chemical and electric propellant on board can be entered instead of a thrust-to-weight and an acceleration, in Python (`Mission(mass0_kg=...)`), MATLAB (`'mass0_kg'`) and the web, behind a Vehicle toggle in both planners; the generic inputs stay the default. The thrust and mass set the thrust-to-weight, the electric thrust is 2ηP/vₑ and stays constant as the mass falls (the spiral time is vₑ/a₀ (1 − e^(−Δv/vₑ))), propellant reads in kg, and an option needing more propellant than is on board is ruled out. The GMAT script carries the spacecraft and a chemical tank. New example: Small tank.
- **Orbital elements:** each orbit's inclination, RAAN and argument of periapsis can be entered instead of an inclination change, in Python (`Mission(elements=...)`), MATLAB (`'elements'`) and the web, behind an Orientation toggle in both planners; the inclination change stays the default. The models get the angle between the planes and each periapsis measured from the node line; apse-to-apse and Lambert accept a start periapsis off the node line. The GMAT script and JSON hold the orbits as entered, in EarthMJ2000Eq, and the script also reports RAAN and argument of periapsis. New example: Kourou GTO (i 6°, ω 178°) to GEO.
- **GMAT export:** the web planner and the MATLAB app save the selected transfer as a NASA GMAT script and a JSON file of its burns; Python has `gmat_script` and `burns_json`. The script flies the burn plan from the start state with the tool's physics, including its Earth μ, and reports the orbit it arrives on.
- **Exact burn plans for every impulsive transfer:** the circular transfers now also carry a start state and inertial Δv vectors, built with an analytic Kepler propagator.
- **Verification:** GMAT flies the exported scripts and arrives on target to 1.1 × 10⁻⁹; the export and the orbital elements are tested in all three engines.

### Changed

- **Circular transfers' burn plans** now turn the plane toward the target plane itself and change plane on the node line, so they arrive in the right plane, not only at the right inclination. A raising one-tangent transfer starts earlier on the circular orbit for this. Δv and flight times are unchanged.
- **Kepler propagator:** hyperbolic arcs start Newton from Vallado's guess; the elliptic one overflowed on a fast 59 km/s Lambert arc.
- **Test suites:** `validate.py` 108 checks, `odst_tests.m` 126, `test_engine.js` 121; the GMAT cross-check runs 36 checks, 33 of them on exported scripts.
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
