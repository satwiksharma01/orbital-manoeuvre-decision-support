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
// ---- elliptical orbits and Lambert: the checks of python/validate.py section 5
const check = (name, got, ref, tol) => {
  const ok = Math.abs(got - ref) <= tol; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] ${name}: got ${got.toFixed(5)}, ref ${ref.toFixed(5)} (tol ${tol})`);
};
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
  [{ alt1: 250, apo1: 35786, alt2: 35786, di_deg: 28.5 }, [[1.833344084, 0], [1.833344097, 20504.3295]], "apse"],
  [{ alt1: 400, alt2: 600, apo2: 39700, dw_deg: 270, di_deg: 0 }, [[2.476172306, 21378.6762], [2.476173769, 21646.7968]], "apse"],
  [{ alt1: 500, apo1: 5000, alt2: 800, apo2: 12000, dw_deg: 60, di_deg: 0 }, [null, [1.183670980, 9359.7970]], "lambert"],
  [{ alt1: 500, apo1: 5000, alt2: 800, apo2: 12000, dw_deg: 60, di_deg: 10 }, [null, [1.859013203, 3601.2405]], "lambert"],
];
for (const [x, refs, rec] of ell) {
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
    ok = Math.abs(o.dv - dv) < 1e-6 && Math.abs(o.tof - tof) < 1; fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] ${tag}: ${o.key} dv ${o.dv.toFixed(9)} vs Python ${dv.toFixed(9)} ` +
                `(${Math.abs(o.dv - dv).toExponential(1)} km/s), TOF ${o.tof.toFixed(4)} vs ${tof.toFixed(4)} s`);
    const sim = O.simulateVec(o.y0, o.plan, 3600, 20000);
    const [dh, de] = O.targetMismatch(sim.yf, r.r2, r.ra2, m.di_deg * Math.PI / 180, m.dw_deg * Math.PI / 180);
    ok = dh < 1e-6 && de < 1e-6 && Math.abs(sim.dvUsed - o.dv) < 1e-9; fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] ${tag}: ${o.key} propagated (RK4) arrival |dh|/h = ${dh.toExponential(1)}, |de| = ${de.toExponential(1)}`);
  });
}
console.log(fails ? `${fails} FAILED` : "ALL CHECKS PASSED");
process.exit(fails ? 1 : 0);
