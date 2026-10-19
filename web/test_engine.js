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
console.log(fails ? `${fails} FAILED` : "ALL CHECKS PASSED");
process.exit(fails ? 1 : 0);
