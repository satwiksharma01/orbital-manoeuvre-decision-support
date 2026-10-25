// Extracts the engine from index.html and checks it against python/odst.py reference values.
// Run: node web/test_engine.js
const fs = require("fs"), path = require("path");
const html = fs.readFileSync(path.join(__dirname, "index.html"), "utf8");
const src = html.match(/<script id="engine">([\s\S]*?)<\/script>/)[1];
const module_ = { exports: {} };
new Function("module", src)(module_);
const O = module_.exports;
const base = { alt1: 400, alt2: 35786, di_deg: 28.5, dv_budget: 6.5, margin_pct: 5, tof_max_h: 4800, propulsion: "both", isp_chem: 320,
  isp_ep: 1600, accel_ep: 5e-4, rb_max: 384400, k_one_tangent: 1.5, weights: O.PRESETS["Balanced"] };
// reference values printed by python/validate.py and odst.py
const expect = [
  [base, { hohmann: 4.22154, hohmann_split: 4.19630, one_tangent: 5.47051, low_thrust: 5.89752 }, "hohmann"],
  [{ ...base, alt1: 7000 - O.RE, alt2: 105000 - O.RE, di_deg: 0, tof_max_h: 1440, propulsion: "chemical", weights: O.PRESETS["Propellant-critical"] },
   { hohmann: 4.04633, bielliptic: 3.99791, one_tangent: 4.85310 }, "hohmann"],
  [{ ...base, alt2: 800, di_deg: 0, dv_budget: 1, tof_max_h: 0.6, weights: O.PRESETS["Time-critical"] }, { one_tangent: 0.32585 }, "one_tangent"],
  [{ ...base, alt1: 35786, alt2: 400, di_deg: 10 }, { hohmann: 3.90496, hohmann_split: 3.90030, one_tangent: 7.26043 }, "hohmann"],
];
let fails = 0;
for (const [m, dvs, rec] of expect) {
  const r = O.evaluate(m);
  for (const [k, v] of Object.entries(dvs)) {
    const got = r.opts.find((o) => o.key === k).dv;
    const ok = Math.abs(got - v) < 1e-4; fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] ${m.alt1.toFixed(0)}->${m.alt2.toFixed(0)} ${k}: ${got.toFixed(5)} (ref ${v})`);
  }
  const ok = r.rec && r.rec.key === rec; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] recommendation ${r.rec && r.rec.key} (ref ${rec})`);
}
for (const [budget, margin, want] of [[6, 0, true], [6, 5, false], [6.5, 5, true]]) {
  const lt = O.evaluate({ ...base, dv_budget: budget, margin_pct: margin }).opts.find((o) => o.key === "low_thrust");
  const ok = lt.feasible === want; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] budget ${budget} km/s, margin ${margin}%: low thrust feasible = ${lt.feasible}`);
}
for (const [tw, ref] of [[0.5, 14.95], [0.1, 373.81]]) {   // python/odst.py, Hohmann LEO -> GEO
  const loss = 1e3 * O.evaluate({ ...base, tw_chem: tw }).opts.find((o) => o.key === "hohmann").dv_loss;
  const ok = Math.abs(loss - ref) < 0.1; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] Hohmann gravity loss at T/W ${tw}: ${loss.toFixed(2)} m/s (ref ${ref})`);
}
const check = (name, got, ref, tol) => {
  const ok = Math.abs(got - ref) <= tol; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] ${name}: got ${got.toFixed(5)}, ref ${ref.toFixed(5)} (tol ${tol})`);
};
{
  const h = O.evaluate({ ...base, tw_chem: 0.2 }).opts.find((o) => o.key === "hohmann");
  check("Propellant fraction includes the gravity loss (Hohmann, T/W 0.2)", h.prop_frac, 1 - Math.exp(-(h.dv + h.dv_loss) / (h.isp * O.G0)), 1e-12);
  check("Same propellant fraction as python/odst.py (Hohmann, T/W 0.2)", h.prop_frac, 0.74717, 1e-5);
  const lt = O.evaluate({ ...base, di_deg: 120 }).opts.find((o) => o.key === "low_thrust");
  fails += lt.applicable;
  console.log(`[${lt.applicable ? "FAIL" : "PASS"}] Edelbaum not offered beyond 114.6 deg (di 120): applicable = ${lt.applicable} (${lt.note})`);
}

// ---- elliptical orbits and Lambert: the checks of python/validate.py section 5
const r1v = [5000, 10000, 2100], r2v = [-14600, 2500, 7000];   // Curtis, Example 5.2
const [v1, v2] = O.lambert(r1v, r2v, 3600);
[-5.9925, 1.9254, 3.2456, -3.3125, -4.1966, -0.38529].forEach((ref, k) =>
  check(`Lambert, Curtis Ex. 5.2, v${k < 3 ? 1 : 2}[${k % 3}] [km/s]`, [...v1, ...v2][k], ref, 1e-4));
const end = O.simulateVec([...r1v, ...v1], [], 3600, 20000).yf;
check("Lambert v1 propagated for 1 h reaches r2 [km]", Math.hypot(...end.slice(0, 3).map((x, i) => x - r2v[i])), 0, 1e-6);
const di20 = 20 * Math.PI / 180;
check("Apse-to-apse between circular orbits == Hohmann, optimal split [km/s]", O.apseTransfer(7000, 7000, 14000, 14000, di20, 0).dv,
      O.hohmann(7000, 14000, di20, O.hohmannSplitOpt(7000, 14000, di20)).dv, 1e-9);
check("Lambert optimum, coplanar 7000 -> 14000 km, vs Hohmann [m/s]", 1e3 * O.lambertOptimal(7000, 7000, 14000, 14000, 0, 0, 1e9).dv,
      1e3 * O.hohmann(7000, 14000, 0).dv, 1.0);
// reference dv [km/s] and TOF [s] from python/odst.py; null = not applicable
const ell = [
  [{ alt1: 250, apo1: 35786, alt2: 35786, di_deg: 28.5 }, [[1.833344084, 0], [1.833344097, 20504.3564]], "apse"],
  [{ alt1: 400, alt2: 600, apo2: 39700, dw_deg: 270, di_deg: 0 }, [[2.476172306, 21378.6762], [2.476173769, 21646.7968]], "apse"],
  [{ alt1: 500, apo1: 5000, alt2: 800, apo2: 12000, dw_deg: 60, di_deg: 0 }, [null, [1.183670980, 9359.7970]], "lambert"],
  [{ alt1: 500, apo1: 5000, alt2: 800, apo2: 12000, dw_deg: 60, di_deg: 10 }, [null, [1.859013203, 3601.2428]], "lambert"],
  // 100 km floor active; on it the optimum's time is less sharply determined, so TOF gets 5 s
  [{ alt1: 300, apo1: 12500, alt2: 350, apo2: 8600, dw_deg: 90, di_deg: 0, tof_max_h: 1950 / 3600 }, [null, [2.546787403, 1629.7306]], "lambert", 5],
];
for (const [x, refs, rec, tolT = 1] of ell) {
  const m = { ...base, dw_deg: 0, ...x }, r = O.evaluate(m);
  const tag = `${m.alt1}x${m.apo1 ?? m.alt1} -> ${m.alt2}x${m.apo2 ?? m.alt2} km di=${m.di_deg} dw=${m.dw_deg}`;
  let ok = r.rec && r.rec.key === rec; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] ${tag}: recommends ${r.rec && r.rec.key} (Python: ${rec})`);
  r.opts.forEach((o, i) => {
    if (!o.applicable) {
      ok = refs[i] === null; fails += !ok;
      console.log(`[${ok ? "PASS" : "FAIL"}] ${tag}: ${o.key} not applicable (${o.note})`);
      return;
    }
    const [dv, tof] = refs[i];
    ok = Math.abs(o.dv - dv) < 1e-6 && Math.abs(o.tof - tof) < tolT; fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] ${tag}: ${o.key} dv ${o.dv.toFixed(9)} vs Python ${dv.toFixed(9)} ` +
                `(${Math.abs(o.dv - dv).toExponential(1)} km/s), TOF ${o.tof.toFixed(4)} vs ${tof.toFixed(4)} s`);
    const sim = O.simulateVec(o.y0, o.plan, 3600, 20000);
    const [dh, de] = O.targetMismatch(sim.yf, r.r2, r.ra2, m.di_deg * Math.PI / 180, m.dw_deg * Math.PI / 180);
    ok = dh < 1e-6 && de < 1e-6 && Math.abs(sim.dvUsed - o.dv) < 1e-9; fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] ${tag}: ${o.key} propagated (RK4) arrival |dh|/h = ${dh.toExponential(1)}, |de| = ${de.toExponential(1)}`);
    if (o.key === "lambert") {                  // the arc must stay above the 100 km floor
      const v1 = o.y0.slice(3).map((v, i) => v + o.plan[0][1][i]), r1 = o.y0.slice(0, 3);
      const low = O.arcMinRadius(r1, v1, O.simulateVec([...r1, ...v1], [], o.plan[1][0], 20000).yf.slice(0, 3)) - O.RE;
      ok = low >= O.H_ARC_MIN - 1e-6; fails += !ok;
      console.log(`[${ok ? "PASS" : "FAIL"}] ${tag}: lambert arc lowest altitude ${low.toFixed(4)} km (floor ${O.H_ARC_MIN} km)`);
    }
  });
}
// ---- export: exact burn plans, GMAT script and JSON (as python/validate.py section 5b)
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const nrm = (a) => Math.hypot(...a);
let worst = 0;      // Kepler propagation vs Kepler's equation, [rp, ra, E0, E1]
for (const [rp, ra, E0, E1] of [[6778, 6778, 0.3, 2.0], [7000, 42164, 1.0, 4.5], [6678, 384400, 5.0, 7.0], [8000, 30000, 2.5, 3.0]]) {
  const a = (rp + ra) / 2, e = (ra - rp) / (ra + rp);
  const th = (E) => 2 * Math.atan2(Math.sqrt(1 + e) * Math.sin(E / 2), Math.sqrt(1 - e) * Math.cos(E / 2));
  const t = (E) => Math.sqrt(a ** 3 / O.MU) * (E - e * Math.sin(E));
  const s0 = O.orbitState(rp, ra, 0.4, 1.1, th(E0)), s1 = O.orbitState(rp, ra, 0.4, 1.1, th(E1));
  const y = O.kepler([...s0.r, ...s0.v], t(E1) - t(E0));
  worst = Math.max(worst, nrm(y.slice(0, 3).map((x, i) => x - s1.r[i])) / nrm(s1.r), nrm(y.slice(3).map((x, i) => x - s1.v[i])) / nrm(s1.v));
}
check("Kepler propagation vs Kepler's equation, 4 arcs up to e = 0.97 [1e-12 relative]", 1e12 * worst, 0, 100);
const circ = [...expect.map((x) => x[0]), { ...base, alt2: 400, di_deg: 60 }, { ...base, alt1: 35786, alt2: 400, tof_max_h: 3 }];
for (const m of circ) {
  const r = O.evaluate(m);
  for (const o of r.opts) {
    if (!o.cplan || !o.applicable) continue;
    const sim = O.simulateVec(o.y0, o.plan, 3600, 200000), y = sim.yf, rv = y.slice(0, 3), v = y.slice(3), h = cross(rv, v);   // fine steps: coasts out to e = 0.97
    const a = 1 / (2 / nrm(rv) - (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) / O.MU);
    const e = nrm(cross(v, h).map((x, i) => x / O.MU - rv[i] / nrm(rv))), inc = Math.acos(h[2] / nrm(h)) * 180 / Math.PI;
    const ok = Math.abs(a - r.r2) / r.r2 < 1e-8 && e < 1e-6 && Math.abs(inc - m.di_deg) < 1e-6 && Math.abs(sim.dvUsed - o.dv) < 1e-12; fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] ${m.alt1.toFixed(0)}->${m.alt2.toFixed(0)} km di=${m.di_deg}: ${o.key} vector plan (RK4) ` +
                `da/a=${((a - r.r2) / r.r2).toExponential(1)} e=${e.toExponential(1)} di=${(inc - m.di_deg).toExponential(1)} deg, sum|dv| - dv=${(sim.dvUsed - o.dv).toExponential(1)}`);
  }
}
// same plans as python/odst.py: closed-form burns (no optimiser), [t, [dvx, dvy, dvz]]
const pyRefV = { hohmann: [[0, [0, 2.3974725168152418, 0]], [19048.562509797113, [9.211582468598876e-16, -1.0838903670048219, 1.467101568882389]]],
                 one_tangent: [[0, [0, 2.638184468665157, 0]], [10862.249834254219, [1.426823308614846, -1.958020944631118, 1.4671015688823885]]] };
for (const [k, ref] of Object.entries(pyRefV)) {
  const o = O.evaluate(base).opts.find((x) => x.key === k);
  const d = Math.max(...ref.flatMap(([t, dv], i) => [t, ...dv].map((x, j) => Math.abs([o.plan[i][0], ...o.plan[i][1]][j] - x) / Math.max(1, Math.abs(x)))));
  const ok = d < 1e-12; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] LEO -> GEO ${k} vector plan vs Python: largest difference ${d.toExponential(1)}`);
}
// GMAT script and JSON read back to the exact plan (GMAT flies them: python/gmat_crosscheck.py)
for (const m of [base, { ...base, alt1: 400, alt2: 600, apo2: 39700, dw_deg: 270, di_deg: 0 }]) {
  for (const o of O.evaluate(m).opts) {
    if (!o.plan || !o.plan.length) continue;
    const txt = O.gmatScript(m, o), js = JSON.parse(O.burnsJson(m, o));
    const val = (k) => +txt.match(new RegExp(`^${k.replace(/\./g, "\\.")} = (\\S+);$`, "m"))[1];
    const plan = []; let t = 0;
    for (const line of txt.slice(txt.indexOf("BeginMissionSequence")).split("\n")) {
      if (line.startsWith("Propagate") && plan.length < o.plan.length) t += +line.match(/ElapsedSecs = (\S+)}/)[1];
      else if (line.startsWith("Maneuver")) { const b = plan.length + 1; plan.push([t, [1, 2, 3].map((j) => val(`Burn${b}.Element${j}`))]); }
    }
    const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
    const ok = val("Earth.Mu") === O.MU && same(["X", "Y", "Z", "VX", "VY", "VZ"].map((k) => val("Sat." + k)), o.y0.map((x) => x + 0)) &&
               same(js.start_state_km_kms, o.y0.map((x) => x + 0)) && plan.length === o.plan.length &&
               plan.every(([tb, dv], i) => Math.abs(tb - o.plan[i][0]) <= 1e-9 * Math.max(1, o.plan[i][0]) && same(dv, o.plan[i][1].map((x) => x + 0))) &&
               same(js.burns.map((b) => [b.t_s, b.dv_kms]), o.plan.map(([tb, dv]) => [tb, dv.map((x) => x + 0)]));
    fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] ${o.key}: GMAT script and JSON give back the plan (${plan.length} burns, start state, Earth mu)`);
  }
}
console.log(fails ? `${fails} FAILED` : "ALL CHECKS PASSED");
process.exit(fails ? 1 : 0);
