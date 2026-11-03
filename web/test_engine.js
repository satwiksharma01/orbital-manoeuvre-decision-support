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
worst = 0;          // ... and hyperbolic arcs (fast Lambert transfers): M = e sinh F - F; [rp, e, F0, F1]
for (const [rp, e, F0, F1] of [[7000, 1.05, -1, 2], [12000, 3, -2.5, 0.5], [20000, 7.5, 0.2, 2.9]]) {
  const a = rp / (1 - e), p = rp * (1 + e), { P, Q } = O.perifocal(0.4, 1.1, 0.7);
  const th = (F) => 2 * Math.atan(Math.sqrt((e + 1) / (e - 1)) * Math.tanh(F / 2)), t = (F) => Math.sqrt((-a) ** 3 / O.MU) * (e * Math.sinh(F) - F);
  const st = (nu) => [...P.map((c, i) => p / (1 + e * Math.cos(nu)) * (Math.cos(nu) * c + Math.sin(nu) * Q[i])),
                      ...P.map((c, i) => Math.sqrt(O.MU / p) * (-Math.sin(nu) * c + (e + Math.cos(nu)) * Q[i]))];
  const y = O.kepler(st(th(F0)), t(F1) - t(F0)), ref = st(th(F1));
  worst = Math.max(worst, nrm(y.slice(0, 3).map((x, i) => x - ref[i])) / nrm(ref.slice(0, 3)), nrm(y.slice(3).map((x, i) => x - ref[i + 3])) / nrm(ref.slice(3)));
}
check("Kepler propagation vs the hyperbolic Kepler equation, 3 arcs, e 1.05 to 7.5 [1e-12 relative]", 1e12 * worst, 0, 100);
const circ = [...expect.map((x) => x[0]), { ...base, alt2: 400, di_deg: 60 }, { ...base, alt1: 35786, alt2: 400, tof_max_h: 3 }];
for (const m of circ) {
  const r = O.evaluate(m);
  for (const o of r.opts) {
    if (!o.cplan || !o.applicable) continue;
    const sim = O.simulateVec(o.y0, o.plan, 3600, 200000);   // fine steps: coasts out to e = 0.97
    // the target plane itself (normal [0, -sin di, cos di]), not only its inclination
    const [dh, de] = O.targetMismatch(sim.yf, r.r2, r.ra2, m.di_deg * Math.PI / 180, 0);
    const ok = dh < 1e-6 && de < 1e-6 && Math.abs(sim.dvUsed - o.dv) < 1e-12; fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] ${m.alt1.toFixed(0)}->${m.alt2.toFixed(0)} km di=${m.di_deg}: ${o.key} vector plan (RK4) ` +
                `|dh|/h=${dh.toExponential(1)} |de|=${de.toExponential(1)}, sum|dv| - dv=${(sim.dvUsed - o.dv).toExponential(1)}`);
  }
}
// same plans as python/odst.py: closed-form burns (no optimiser), start state and [t, [dvx, dvy, dvz]]
const pyRefV = { hohmann: [[6778.137, 0, 0, 0, 7.668558175407055, 0],
                           [[0, [0, 2.3974725168152418, 0]], [19048.562509797113, [9.211582468598876e-16, -1.0838903670048219, -1.467101568882389]]]],
                 one_tangent: [[-5964.52423775776, -3219.874501897429, 0, 3.6428586402591523, -6.748063827253831, 0],
                               [[0, [1.2532385966758652, -2.3215103511524484, 0]], [10862.249834254219, [-2.18568965499367, 1.0451944822253934, 1.467101568882389]]]] };
const rel = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(b));
for (const [k, [y0, ref]] of Object.entries(pyRefV)) {
  const o = O.evaluate(base).opts.find((x) => x.key === k);
  const d = Math.max(...y0.map((x, j) => rel(o.y0[j], x)),
                     ...ref.flatMap(([t, dv], i) => [t, ...dv].map((x, j) => rel([o.plan[i][0], ...o.plan[i][1]][j], x))));
  const ok = d < 1e-12; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] LEO -> GEO ${k} start state and vector plan vs Python: largest difference ${d.toExponential(1)}`);
}
// GMAT script and JSON read back to the exact plan (GMAT flies them: python/gmat_crosscheck.py)
for (const m of [base, { ...base, alt1: 400, alt2: 600, apo2: 39700, dw_deg: 270, di_deg: 0 },
                 { ...base, alt1: 400, alt2: 600, apo2: 39700, elements: [51.6, 30, 0, 63.4, 100, 270], tof_max_h: 30 }]) {
  for (const o of O.evaluate(m).opts) {
    if (!o.plan || !o.plan.length) continue;
    const txt = O.gmatScript(m, o), js = JSON.parse(O.burnsJson(m, o)), ip = O.inertialPlan(m, o);
    const val = (k) => +txt.match(new RegExp(`^${k.replace(/\./g, "\\.")} = (\\S+);$`, "m"))[1];
    const plan = []; let t = 0;
    for (const line of txt.slice(txt.indexOf("BeginMissionSequence")).split("\n")) {
      if (line.startsWith("Propagate") && plan.length < ip.plan.length) t += +line.match(/ElapsedSecs = (\S+)}/)[1];
      else if (line.startsWith("Maneuver")) { const b = plan.length + 1; plan.push([t, [1, 2, 3].map((j) => val(`Burn${b}.Element${j}`))]); }
    }
    const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
    const ok = val("Earth.Mu") === O.MU && same(["X", "Y", "Z", "VX", "VY", "VZ"].map((k) => val("Sat." + k)), ip.y0.map((x) => x + 0)) &&
               same(js.start_state_km_kms, ip.y0.map((x) => x + 0)) && plan.length === ip.plan.length &&
               plan.every(([tb, dv], i) => Math.abs(tb - ip.plan[i][0]) <= 1e-9 * Math.max(1, ip.plan[i][0]) && same(dv, ip.plan[i][1].map((x) => x + 0))) &&
               same(js.burns.map((b) => [b.t_s, b.dv_kms]), ip.plan.map(([tb, dv]) => [tb, dv.map((x) => x + 0)])) &&
               same(js.elements_deg, m.elements ?? null);
    fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] ${o.key}: GMAT script and JSON give back the plan (${plan.length} burns, start state, Earth mu${m.elements ? ", elements" : ""})`);
  }
}

// ---- orbital elements (as python/validate.py section 5c)
const dvs = (m) => Object.fromEntries(O.evaluate(m).opts.filter((o) => o.applicable).map((o) => [o.key, o.dv]));
const diff = (a, b) => (same2(Object.keys(a), Object.keys(b)) ? Math.max(...Object.keys(a).map((k) => Math.abs(a[k] - b[k]))) : Infinity);
const same2 = (a, b) => JSON.stringify(a) === JSON.stringify(b);
for (const [relm, el] of [[{ ...base, alt1: 400, alt2: 600, apo2: 39700, di_deg: 20, dw_deg: 270 }, [0, 0, 0, 20, 0, 270]],
                          [{ ...base, alt1: 250, apo1: 35786, alt2: 35786, di_deg: 28.5 }, [0, 0, 0, 28.5, 0, 0]]])
  check(`Elements [${el}] == relative inputs, every option [km/s]`, diff(dvs({ ...relm, elements: el }), dvs(relm)), 0, 1e-12);
const ell2 = { ...base, alt1: 300, apo1: 12000, alt2: 800, apo2: 20000, tof_max_h: 20 };
check("Both orbits turned 40 deg about the pole: Lambert, 3D [km/s]",
      diff(dvs({ ...ell2, elements: [28.5, 10, 30, 40, 70, 120] }), dvs({ ...ell2, elements: [28.5, 50, 30, 40, 110, 120] })), 0, 1e-9);
const raan = O.resolve({ ...base, alt1: 500, alt2: 500, elements: [51.6, 0, 0, 51.6, 10, 0] });
const refDi = 2 * Math.asin(Math.sin(51.6 * Math.PI / 180) * Math.sin(5 * Math.PI / 180));
check("RAAN 0 -> 10 deg at i = 51.6 deg: angle between the planes [deg]", raan.di_deg, refDi * 180 / Math.PI, 1e-12);
check("  its single-burn plane change, 2 v sin(di/2) [km/s]", dvs(raan).hohmann, 2 * Math.sqrt(O.MU / (O.RE + 500)) * Math.sin(refDi / 2), 1e-12);
const gto = { ...base, alt1: 250, apo1: 35786, alt2: 35786 };
const g180 = dvs({ ...gto, elements: [6, 0, 180, 0, 0, 0] });
check("GTO i 6, argp 180 -> GEO: apse-to-apse == relative di 6 [km/s]", g180.apse, dvs({ ...gto, di_deg: 6 }).apse, 1e-12);
{
  const o = O.evaluate({ ...gto, elements: [6, 0, 178, 0, 0, 0] }).opts;
  const ok = !o[0].applicable && o[1].applicable && o[1].dv > g180.apse; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] GTO i 6, argp 178 -> GEO: apse-to-apse not offered (${o[0].note}); Lambert ${o[1].dv.toFixed(6)} km/s, ` +
              `${(1e3 * (o[1].dv - g180.apse)).toFixed(1)} m/s above the aligned case`);
}
// same derived geometry and results as python/odst.py: [di, w1, dw] deg, then [dv km/s, TOF s] per applicable option
for (const [x, geo, refs] of [
  [{ alt1: 400, alt2: 35786, elements: [51.6, 120, 0, 0, 0, 0] }, [51.60000000000001, 0, 0],
   [[4.824649844, 19048.5625], [4.784537299, 19048.5625], [4.501603228, 920556.7324], [5.901213052, 10862.2498], [7.805548982, 15611097.9632]]],
  [{ alt1: 250, apo1: 35786, alt2: 35786, elements: [6, 0, 178, 0, 0, 0] }, [6.000000000000001, 358, 0], [[1.496760716, 3479.8681]]],
  [{ alt1: 400, alt2: 600, apo2: 39700, elements: [51.6, 30, 0, 63.4, 100, 270], tof_max_h: 30 }, [58.815692017167976, 0, 210.591986530302], [[4.442609368, 8948.9431]]]]) {
  const r = O.evaluate({ ...base, ...x }), opts = r.opts.filter((o) => o.applicable);
  const ok = Math.max(...[r.m.di_deg, r.m.w1_deg, r.m.dw_deg].map((v, i) => Math.abs(v - geo[i]))) < 1e-9 && opts.length === refs.length &&
             opts.every((o, i) => Math.abs(o.dv - refs[i][0]) < 1e-6 && Math.abs(o.tof - refs[i][1]) < 1); fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] [${x.elements}]: di ${r.m.di_deg.toFixed(6)}, w1 ${r.m.w1_deg.toFixed(6)}, dw ${r.m.dw_deg.toFixed(6)} deg ` +
              `and ${opts.length} options as Python (dv to 1e-6 km/s, TOF to 1 s)`);
}
// every plan, turned into EarthMJ2000Eq and propagated, arrives on the target given by its own elements
for (const x of [{ alt1: 400, alt2: 35786, elements: [51.6, 120, 0, 0, 0, 0] }, { alt1: 400, alt2: 35786, elements: [28.5, 0, 0, 0, 0, 0], tof_max_h: 3 },
                 { alt1: 500, alt2: 500, elements: [51.6, 0, 0, 51.6, 10, 0] }, { alt1: 250, apo1: 35786, alt2: 35786, elements: [6, 0, 178, 0, 0, 0] },
                 { alt1: 400, alt2: 600, apo2: 39700, elements: [51.6, 30, 0, 63.4, 100, 270], tof_max_h: 30 },
                 { alt1: 300, apo1: 12000, alt2: 800, apo2: 20000, elements: [28.5, 10, 30, 40, 70, 120], tof_max_h: 20 }]) {
  const m = { ...base, ...x }, r = O.evaluate(m), { P, W } = O.perifocal(...m.elements.slice(3).map((v) => v * Math.PI / 180));
  const rp = r.r2, ra = r.ra2, p = 2 * rp * ra / (rp + ra), hn = Math.sqrt(O.MU * p);
  for (const o of r.opts) {
    if (!o.plan || !o.plan.length || !o.applicable) continue;
    const ip = O.inertialPlan(m, o), sim = O.simulateVec(ip.y0, ip.plan, 3600, 200000), rv = sim.yf.slice(0, 3), v = sim.yf.slice(3), h = cross(rv, v);
    const dh = nrm(h.map((c, i) => c - hn * W[i])) / hn;
    const de = nrm(cross(v, h).map((c, i) => c / O.MU - rv[i] / nrm(rv) - (ra - rp) / (ra + rp) * P[i]));
    const ok = dh < 1e-6 && de < 1e-6 && Math.abs(sim.dvUsed - o.dv) < 1e-9; fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] [${m.elements}] ${o.key}: dv ${o.dv.toFixed(6)}, arrival (RK4) vs the target's elements ` +
                `|dh|/h=${dh.toExponential(1)} |de|=${de.toExponential(1)}`);
  }
}
// ---- spacecraft model (as python/validate.py section 4c)
const scBase = { ...base, thrust_chem_N: 4900, prop_chem_kg: 800, power_ep_kW: 6.5, eff_ep: 0.6, prop_ep_kg: 100 };
{
  const G = O.G0 * 1e3, gen = O.evaluate(base).opts, scm = { ...scBase, mass0_kg: 1000, thrust_chem_N: 0.5 * 1000 * G }, sc = O.evaluate(scm).opts;
  const d = Math.max(...gen.map((a, i) => (a.key === "low_thrust" ? 0 : Math.abs(a.dv - sc[i].dv) + Math.abs(a.dv_loss - sc[i].dv_loss) + Math.abs(a.tof - sc[i].tof))));
  check("Spacecraft with thrust = 0.5 m0 g0 vs generic T/W 0.5: chemical options [dv + loss + TOF difference]", d, 0, 1e-9);
  check("  propellant = m0 x propellant fraction, Hohmann [kg]", sc[0].prop_kg, 1000 * gen[0].prop_frac, 1e-9);
  const F = O.epThrust(scm);
  check("Electric thrust 2 eta P / ve, 6.5 kW, eta 0.6, Isp 1600 s [N]", F, 2 * 0.6 * 6500 / (1600 * G), 1e-12);
  // Edelbaum's averaged equations in time with a = F/m(t) and the yaw law in the velocity increment u (RK4, fine steps)
  const lt = sc.find((o) => o.key === "low_thrust"), r1 = O.RE + 400, r2 = O.RE + 35786, di = 28.5 * Math.PI / 180;
  const V0 = Math.sqrt(O.MU / r1), V1 = Math.sqrt(O.MU / r2), b0 = Math.atan2(Math.sin(Math.PI / 2 * di), V0 / V1 - Math.cos(Math.PI / 2 * di));
  const ve = 1600 * O.G0, f0 = F * 1e-3 / 1000;
  const f = ([V, , mr, u]) => { const b = Math.atan2(V0 * Math.sin(b0), V0 * Math.cos(b0) - u); return [-f0 / mr * Math.cos(b), 2 / Math.PI * f0 / mr * Math.sin(b) / V, -f0 / ve, f0 / mr]; };
  let y = [V0, 0, 1, 0]; const n = 200000, h = lt.tof / n, ax = (a, k, c) => a.map((v, i) => v + k * c[i]);
  for (let i = 0; i < n; i++) { const k1 = f(y), k2 = f(ax(y, h / 2, k1)), k3 = f(ax(y, h / 2, k2)), k4 = f(ax(y, h, k3)); y = y.map((v, j) => v + h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])); }
  check("Constant-thrust spiral LEO -> GEO, di 28.5: speed after the closed-form TOF [km/s]", y[0], V1, 1e-8);
  check("  plane change after the closed-form TOF [deg]", y[1] * 180 / Math.PI, 28.5, 1e-6);
  check("  mass left, m/m0, vs exp(-dv/ve)", y[2], Math.exp(-lt.dv / ve), 1e-9);
  check("Constant thrust with Isp -> infinity becomes constant acceleration [TOF ratio]", O.lowThrust(r1, r2, di, 5e-4, 1e9).tof / O.lowThrust(r1, r2, di, 5e-4).tof, 1, 1e-6);
  // the propellant on board is a hard limit, at the same dv as the budget check
  const need = 1000 * (1 - Math.exp(-(sc[0].dv + sc[0].dv_loss) * 1.05 / (sc[0].isp * O.G0)));
  let ok = O.evaluate({ ...scm, prop_chem_kg: need + 1e-6 }).opts[0].feasible && !O.evaluate({ ...scm, prop_chem_kg: need - 1e-6 }).opts[0].feasible; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] Hohmann needs ${need.toFixed(6)} kg with the 5% margin: feasible with 1e-6 kg more on board, not with 1e-6 kg less`);
  const t = O.evaluate({ ...scBase, mass0_kg: 1000, prop_chem_kg: 600, prop_ep_kg: 350 });
  ok = t.rec && t.rec.key === "low_thrust" && t.opts.every((o) => o.key === "low_thrust" || !o.feasible); fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] 600 kg of chemical propellant on board: every chemical option is out (${t.opts[0].reasons[0]}); ` +
              `low thrust is recommended, ${t.rec.prop_kg.toFixed(0)} kg in ${(t.rec.tof / 86400).toFixed(1)} days`);
  // same as python/odst.py, default spacecraft, LEO -> GEO: [dv km/s, loss km/s, TOF s, propellant kg]
  const ref = [[4.221535383, 0.014972551, 19048.5625, 740.762642], [4.19629709, 0.015327503, 19048.5625, 738.698885],
               [5.470508212, 0.01867395, 10862.2498, 826.085245], [5.897522253, 0, 9888981.43, 313.303212]];
  const got = O.evaluate({ ...scBase, mass0_kg: 1000 }).opts.filter((o) => o.applicable).map((o) => [o.dv, o.dv_loss, o.tof, o.prop_kg]);
  ok = got.length === ref.length && got.every((g, i) => [0, 1, 3].every((j) => Math.abs(g[j] - ref[i][j]) < 1e-6) && Math.abs(g[2] - ref[i][2]) < 1e-3); fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] Default spacecraft, LEO -> GEO: dv, loss, TOF and propellant of ${got.length} options as Python`);
  // the export carries the spacecraft
  const txt = O.gmatScript(scm, sc[0]), js = JSON.parse(O.burnsJson(scm, sc[0]));
  ok = txt.includes("Sat.DryMass = 200;") && txt.includes("Fuel.FuelMass = 800;") && txt.split("DecrementMass = true;").length === 3 &&
       js.spacecraft.mass0_kg === 1000 && js.spacecraft.propellant_kg === sc[0].prop_kg; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] GMAT script with the spacecraft's dry mass, an 800 kg tank and mass-decrementing burns; JSON with the propellant, ${sc[0].prop_kg.toFixed(1)} kg`);
}
// ---- uncertainty: Monte Carlo over Isp, thrust and pointing errors; the Pareto front (as python/validate.py section 4d)
{
  const u = O.mrg32k3a(5, [0, 0, 1234567, 0, 0, 1234567]).map((x) => Math.round(x / 2.328306549295727688e-10) + 1);
  let ok = [1459213977, 2827710106, 4245671317, 3877608661, 2595287583].every((v, i) => u[i] === v); fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] MRG32k3a reproduces the published test vector (seed 1234567, first five outputs)`);
  // the same normals as Python: first and last run
  const Z = O.mcNormals(), py = [[-0.847924823347079, 1.8460727873862615, 0.7028567229701457, -1.361475967116543, -1.6978660974898043, -0.3615845453046289, -0.7451327634504835, 0.953673039659349],
    [0.3355242415679864, 0.630794942430949, -1.095285273790694, 1.0227744211709222, 0.7098950824471811, 0.858224862332908, -0.7171417966521667, -0.10778190475429406]];
  check(`${O.MC_RUNS} runs: the first 8 normals of the first and last run vs Python [largest difference]`,
        Math.max(...[0, 1].flatMap((j) => Z[j ? O.MC_RUNS - 1 : 0].slice(0, 8).map((z, i) => Math.abs(z - py[j][i])))), 0, 1e-12);
  check("  sum of the first 8 normals of every run (32000) vs Python", Z.flatMap((z) => z.slice(0, 8)).reduce((a, b) => a + b, 0), -224.82505288590977, 1e-9);
  // with every error zero, each run is the nominal one
  const zero = { uncertainty: true, sigma_isp_pct: 0, sigma_thrust_pct: 0, sigma_point_deg: 0, confidence_pct: 90 };
  ok = true;
  for (const m of [base, { ...scBase, mass0_kg: 1000, prop_chem_kg: 500, prop_ep_kg: 350 }, { ...base, alt1: 250, apo1: 35786, propulsion: "chemical" },
                   { ...base, weights: [0.25, 0.55, 0.05, 0.15], tof_max_h: 3300 }]) {
    const a = O.evaluate(m).opts, b = O.evaluate({ ...m, ...zero }).opts;
    ok = ok && a.every((o, i) => o.feasible === b[i].feasible) && b.every((o) => !o.feasible || o.chance === 1);
  }
  fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] Errors set to zero, 4 missions: every feasible option fits in all runs, and feasibility is unchanged`);
  // the same chances as Python, run for run (runs that fit, of 4000)
  const err = { uncertainty: true, sigma_isp_pct: 1, sigma_thrust_pct: 2, sigma_point_deg: 0.5, confidence_pct: 90 };
  const cases = [
    ["Spacecraft near its tank limits", { ...scBase, ...err, mass0_kg: 1000, prop_chem_kg: 765, prop_ep_kg: 330, tof_max_h: 2800, confidence_pct: 0 },
     { hohmann: 3328, hohmann_split: 3740, low_thrust: 3078 }],
    ["Propellant-critical, 3300 h", { ...base, ...err, weights: [0.25, 0.55, 0.05, 0.15], tof_max_h: 3300 }, { hohmann: 4000, hohmann_split: 4000, one_tangent: 4000, low_thrust: 2577 }],
    ["Small tank, 2800 h", { ...scBase, ...err, mass0_kg: 1000, prop_chem_kg: 500, prop_ep_kg: 350, tof_max_h: 2800 }, { low_thrust: 3325 }],
    ["GTO -> GEO, budget 1.95 km/s", { ...base, ...err, alt1: 250, apo1: 35786, propulsion: "chemical", dv_budget: 1.95 }, { apse: 2725, lambert: 2708 }],
    ["LEO -> GEO, budget 4.5 km/s, pointing 1.5 deg", { ...base, ...err, dv_budget: 4.5, sigma_point_deg: 1.5 }, { hohmann: 114, hohmann_split: 440 }],
  ];
  for (const [name, m, ref] of cases) {
    const got = Object.fromEntries(O.evaluate(m).opts.filter((o) => !isNaN(o.chance)).map((o) => [o.key, Math.round(o.chance * O.MC_RUNS)]));
    ok = JSON.stringify(got) === JSON.stringify(ref); fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] ${name}: runs that fit ${JSON.stringify(got)} (Python ${JSON.stringify(ref)})`);
  }
  // the rule: the spiral makes the 3300 h limit in 64 % of runs, so it is out and Hohmann is recommended
  const a = O.evaluate({ ...base, weights: [0.25, 0.55, 0.05, 0.15], tof_max_h: 3300 }), b = O.evaluate(cases[1][1]), lt = b.opts.find((o) => o.key === "low_thrust");
  ok = a.rec.key === "low_thrust" && !lt.feasible && b.rec.key === "hohmann" && lt.reasons.at(-1) === "fits the limits in 64% of runs with the errors, 90% required"; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] Propellant-critical, 3300 h: ${a.rec.name} without errors; with them: ${lt.reasons.at(-1)}; recommended: ${b.rec.name}`);
  const front = O.evaluate(base).opts.filter((o) => o.pareto).map((o) => o.key);
  ok = JSON.stringify(front) === JSON.stringify(["hohmann_split", "one_tangent", "low_thrust"]); fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] Pareto front of LEO -> GEO: ${front}`);
}
// ---- J2: the planes' drift, the spiral, the nodal drift (as python/validate.py section 4e), against Python's values
{
  const PC = O.PRESETS["Propellant-critical"], sc = { ...scBase, mass0_kg: 1000, prop_ep_kg: 150 };
  const cases = {
    raan10: { ...base, alt1: 500, alt2: 500, elements: [51.6, 0, 0, 51.6, 10, 0], tof_max_h: 384, weights: PC },
    raan10inc: { ...base, alt1: 500, alt2: 500, elements: [51.6, 0, 0, 53, 10, 0], tof_max_h: 384, weights: PC },
    spiral: { ...base, alt1: 400, alt2: 1200, elements: [51.6, 0, 0, 51.6, 0, 0], tof_max_h: 1e5 },
    spiralsc: { ...sc, alt1: 400, alt2: 1200, elements: [51.6, 0, 0, 51.6, 0, 0], tof_max_h: 1e5 },
    slow: { ...base, alt1: 400, alt2: 1200, elements: [51.6, 0, 0, 51.6, 0, 0], tof_max_h: 1e5, accel_ep: 1e-4 },
    meo: { ...base, alt1: 400, alt2: 20200, elements: [28.5, 0, 0, 55, 40, 0], tof_max_h: 1e5 },
    descent: { ...base, alt1: 20200, alt2: 500, elements: [55, 0, 0, 55, 3, 0], tof_max_h: 1e5 },
    geo: { ...base, alt1: 400, alt2: 35786, elements: [51.6, 120, 0, 0, 0, 0] },
    kourou: { ...base, alt1: 250, apo1: 35786, alt2: 35786, elements: [6, 0, 178, 0, 0, 0], propulsion: "chemical" },
  };
  const ref = {"raan10": {"rec": "nodal_drift", "opts": {"hohmann": [1.0399338192964325, 0.0, true, 0.0, null, null, null], "hohmann_split": [1.0399338192964325, 0.0, false, null, null, null, null], "bielliptic": [1.0399338223227859, 5676.9780327835915, false, null, null, null, null], "low_thrust": [2.4318678440890706, 4863735.688178142, true, 0.8002103036610191, null, null, 2], "nodal_drift": [0.3042482984417413, 1382401.3823999988, true, null, 783.388963230549, 15.9322694368591, null]}}, "raan10inc": {"rec": "nodal_drift", "opts": {"hohmann": [1.066152021082588, 0.0, true, 0.0, null, null, null], "hohmann_split": [1.066152021082588, 0.0, false, null, null, null, null], "bielliptic": [1.0661520240892794, 5676.9780327835915, false, null, null, null, null], "low_thrust": [1.9929670917756679, 3985934.183551336, true, 0.32027102790091955, null, null, 2], "nodal_drift": [0.49060539494954797, 1382401.3823999998, true, null, 858.484512771689, 15.931725164999492, null]}}, "spiral": {"rec": "hohmann", "opts": {"hohmann": [0.4157515124549431, 3026.206774746056, true, 1.5377694229656846e-05, null, null, null], "hohmann_split": [0.4157361347607134, 3026.206774746056, false, null, null, null, null], "bielliptic": [0.4157361383783885, 6308.857407789015, false, null, null, null, null], "one_tangent": [2.4961364736241793, 916.4530230408586, true, 6.531710199375595e-07, null, null, null], "low_thrust": [1.7514522899225888, 3502904.579845178, true, 1.335392849131551, null, null, 3], "nodal_drift": [null, null, false, null, null, null, null]}}, "spiralsc": {"rec": "hohmann", "opts": {"hohmann": [0.4157515124549431, 3026.206774746056, true, 1.5377694229656846e-05, null, null, null], "hohmann_split": [0.4157361347607134, 3026.206774746056, false, null, null, null, null], "bielliptic": [0.4157361383783885, 6308.857407789015, false, null, null, null, null], "one_tangent": [2.4961364736241793, 916.4530230408586, true, 6.531710199375595e-07, null, null, null], "low_thrust": [1.736063744045925, 3306034.2183484803, true, 1.3200043032548872, null, null, 3], "nodal_drift": [null, null, false, null, null, null, null]}}, "slow": {"rec": "hohmann", "opts": {"hohmann": [0.4157515124549431, 3026.206774746056, true, 1.5377694229656846e-05, null, null, null], "hohmann_split": [0.4157361347607134, 3026.206774746056, false, null, null, null, null], "bielliptic": [0.4157361383783885, 6308.857407789015, false, null, null, null, null], "one_tangent": [2.4961364736241793, 916.4530230408586, true, 6.531710199375595e-07, null, null, null], "low_thrust": [null, null, false, null, null, null, null], "nodal_drift": [null, null, false, null, null, null, null]}}, "meo": {"rec": "hohmann", "opts": {"hohmann": [4.406918154306862, 10717.717121050175, true, 0.001485716537419357, null, null, null], "hohmann_split": [4.330759229284638, 10717.717121050175, true, 0.0014160622915824206, null, null, null], "bielliptic": [4.405432438024526, 32278.66198453821, false, null, null, null, null], "one_tangent": [5.820344649743336, 6055.36761364423, true, 0.00029958441251221046, null, null, null], "low_thrust": [14.719383419439804, 29438766.83887961, true, 8.254497599477219, null, null, 3], "nodal_drift": [4.006970399058953, 3937656.8604975925, true, null, 400.0, 45.45068452982109, null]}}, "descent": {"rec": "hohmann", "opts": {"hohmann": [3.376835683771683, 10765.949782789343, true, -0.0019006893363155974, null, null, null], "hohmann_split": [3.3761343495235128, 10765.949782789343, true, -0.0016026078272135535, null, null, null], "bielliptic": [3.3787363743786605, 32326.89464629656, false, null, null, null, null], "one_tangent": [6.742805504326066, 8981.535790465627, true, -0.0009324085572508523, null, null, null], "low_thrust": [null, null, false, null, null, null, null], "nodal_drift": [3.3723859323697707, 60628.84854822148, true, null, 20200.0, 0.5771168838591683, null]}}, "geo": {"rec": "hohmann", "opts": {"hohmann": [4.82464984429949, 19048.562509797113, true, 0.0, null, null, null], "hohmann_split": [4.784537299351257, 19048.562509797113, true, 0.0, null, null, null], "bielliptic": [4.501603227731081, 920556.7324006476, true, 0.0, null, null, null], "one_tangent": [5.901213052479644, 10862.249834254219, true, 0.0, null, null, null], "low_thrust": [7.805548981621778, 15611097.963243555, true, 0.0, null, null, 0], "nodal_drift": [null, null, false, null, null, null, null]}}, "kourou": {"rec": "lambert", "opts": {"apse": [null, null, false, null, null, null, null], "lambert": [1.4967607162136738, 3479.8680687829756, true, 0.0, null, null, null]}}, "raan10inc_burns": [[0.0, [0.09606579681241367, 0.0, 0.0], "time", false], [2950.1640720218165, [0.09485322144623076, 0.0, 0.0], "apoapsis", false], [1379451.218327978, [-0.0970403476500783, 0.17900980832924193, 0.0], "node", true], [1382401.3823999998, [-0.09606579681241278, 0.0, 0.0], "periapsis", false]], "descent_burns": [[0.0, [-1.3909970758888304, 0.08910037861377262, -6.938893903907228e-18], "time", false], [10765.949782789343, [-1.982987864746827, -2.7755575615628914e-17, 6.938893903907228e-18], "periapsis", false]]};
  const close = (a, b, tol) => (a == null || Number.isNaN(a) ? b == null : b != null && Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)));
  for (const [name, m] of Object.entries(cases)) {
    const r = O.evaluate({ ...m, j2: true }), want = ref[name], bad = [];
    for (const [k, [dv, tof, app, j2dv, alt, wait, corr]] of Object.entries(want.opts)) {
      const o = r.opts.find((x) => x.key === k);
      if (!o) { bad.push(`${k} missing`); continue; }
      if (o.applicable !== app) bad.push(`${k} applicable ${o.applicable}`);
      if (!app) continue;
      if (!close(o.dv, dv, 1e-9) || !close(o.tof, tof, 1e-9)) bad.push(`${k} dv ${o.dv} tof ${o.tof}`);
      if (j2dv != null && !close(o.j2_dv, j2dv, 1e-9)) bad.push(`${k} J2 dv ${o.j2_dv}`);
      if (alt != null && (!close(o.params.alt_d, alt, 1e-9) || !close(o.params.wait_days, wait, 1e-9))) bad.push(`${k} parking ${o.params.alt_d} km ${o.params.wait_days} d`);
      if (corr != null && o.params.corrections !== corr) bad.push(`${k} corrections ${o.params.corrections}`);
    }
    if ((r.rec ? r.rec.key : null) !== want.rec) bad.push(`recommended ${r.rec && r.rec.key}`);
    const ok = !bad.length; fails += !ok;
    const nd = r.opts.find((o) => o.key === "nodal_drift"), lt = r.opts.find((o) => o.key === "low_thrust");
    console.log(`[${ok ? "PASS" : "FAIL"}] J2, ${name}: ${r.opts.length} options as Python (${r.rec ? r.rec.key : "none"} recommended` +
                `${nd && nd.applicable ? `; nodal drift ${nd.dv.toFixed(4)} km/s, ${nd.params.wait_days.toFixed(2)} d at ${nd.params.alt_d.toFixed(0)} km` : ""}` +
                `${lt && lt.applicable && lt.params.corrections ? `; spiral ${lt.dv.toFixed(4)} km/s, ${lt.params.corrections} corrections` : ""})${ok ? "" : " " + bad.join("; ")}`);
  }
  // from a polar orbit the parking orbits tie on time: the cheapest is taken (Python: parking where it starts)
  {
    const o = O.evaluate({ ...base, alt1: 452, alt2: 1273, elements: [90, 6, 0, 84.2, 235.2, 0], tof_max_h: 20, j2: true }).opts.find((x) => x.key === "nodal_drift");
    const ok = o.applicable && !o.feasible && close(o.dv, 0.9626656591468357, 1e-9) && close(o.tof, 37182873.32233598, 1e-9) && close(o.params.alt_d, 452, 1e-9);
    fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] J2, polar start, time limit out of reach: nodal drift ${o.dv.toFixed(4)} km/s parking at ${o.params.alt_d.toFixed(0)} km, as Python`);
  }
  // the local burns J2 plans fly in: VNB components, events, the node-signed plane change
  for (const [name, key] of [["raan10inc", "nodal_drift"], ["descent", "hohmann"]]) {
    const m = O.resolve({ ...cases[name], j2: true }), o = O.evaluate(m).opts.find((x) => x.key === key);
    const got = O.localBurns(m, o).burns, want = ref[name + "_burns"];
    const ok = got.length === want.length && got.every(([t, dv, at, sg], i) => Math.abs(t - want[i][0]) <= 1e-6 * Math.max(1, t) && at === want[i][2] &&
               sg === want[i][3] && dv.every((x, j) => Math.abs(x - want[i][1][j]) < 1e-12)); fails += !ok;
    console.log(`[${ok ? "PASS" : "FAIL"}] J2, ${name} ${key}: local burns as Python (${got.map((b) => b[2]).join(", ")})`);
  }
  // the export: a J2 field, VNB burns, the node crossing, the target alongside
  const m = { ...cases.raan10inc, j2: true }, o = O.evaluate(m).opts.find((x) => x.key === "nodal_drift");
  const txt = O.gmatScript(m, o), js = JSON.parse(O.burnsJson(m, o));
  let ok = txt.includes("GravityField.Earth.Degree = 2;") && txt.split("Axes = VNB;").length === 6 && txt.includes("{Sat.EarthMJ2000Eq.Z = 0}") &&
           txt.includes("If Sat.EarthMJ2000Eq.VZ > 0") && txt.includes("Create Spacecraft TargetSat;") &&
           js.burns.map((b) => b.at).join() === "time,apoapsis,node,periapsis"; fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] J2 export of the nodal drift: J2 field, 5 VNB burns, the node crossing, TargetSat; JSON events ${js.burns.map((b) => b.at)}`);
  // J2 needs absolute inclinations: with only an inclination change it changes nothing
  const a = O.evaluate({ ...base, alt2: 1200, di_deg: 4, j2: true }).opts.map((x) => [x.key, x.dv]), b = O.evaluate({ ...base, alt2: 1200, di_deg: 4 }).opts.map((x) => [x.key, x.dv]);
  ok = JSON.stringify(a) === JSON.stringify(b); fails += !ok;
  console.log(`[${ok ? "PASS" : "FAIL"}] J2 on with only an inclination change: the two-body options, unchanged`);
}
console.log(fails ? `${fails} FAILED` : "ALL CHECKS PASSED");
process.exit(fails ? 1 : 0);
