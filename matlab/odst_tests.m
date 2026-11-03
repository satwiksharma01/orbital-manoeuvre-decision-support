function nFail = odst_tests()
%ODST_TESTS Validation suite: textbook benchmarks + numerical propagation of
%   every burn plan. Mirrors python/validate.py. Returns the number of failures.
C = odst_const(); RE = C.RE;
nFail = 0;
fprintf('=== 1. Textbook benchmarks ===\n');
[dv, tof] = odst_hohmann(RE + 300, 42164, 0);
check('Hohmann 300 km -> GEO dv [km/s]', dv, 3.893, 0.005);
check('Hohmann 300 km -> GEO TOF [h]', tof/3600, 5.27, 0.02);
r1 = 7000;
Rs = [11 15 16]; expectBi = [false true true];
for i = 1:3
    dh = odst_hohmann(r1, Rs(i)*r1, 0);
    db = odst_bielliptic(r1, Rs(i)*r1, 1e4*Rs(i)*r1, 0);
    ok = (db < dh) == expectBi(i);
    report(ok, sprintf('R = %g: Hohmann %.4f vs bi-elliptic(rB->inf) %.4f', Rs(i), dh, db));
end
check('Simple plane change 30 deg at 7000 km', odst_hohmann(7000, 7000, 30*pi/180), ...
      2*sqrt(C.MU/7000)*sin(15*pi/180), 1e-9);
check('One-tangent(k=1) == Hohmann dv', odst_onetangent(7000, 20000, 1 + 1e-10, 0), ...
      odst_hohmann(7000, 20000, 0), 1e-4);

fprintf('\n=== 2. Numerical verification of burn plans (two-body propagation) ===\n');
cases = {odst_mission('alt1', 400, 'alt2', 35786, 'di_deg', 28.5), ...
         odst_mission('alt1', 7000 - RE, 'alt2', 105000 - RE, 'di_deg', 0, 'tof_max_h', 24*60), ...
         odst_mission('alt1', 400, 'alt2', 800, 'di_deg', 0, 'tof_max_h', 0.6), ...
         odst_mission('alt1', 35786, 'alt2', 400, 'di_deg', 10)};
for c = 1:numel(cases)
    m = cases{c};
    opts = odst_evaluate(m);
    r2 = RE + m.alt2;
    for i = 1:numel(opts)
        o = opts(i);
        if isempty(o.plan) || ~o.applicable, continue; end
        [~, ~, yf, dvSim] = odst_simulate(RE + m.alt1, o.plan, 3600);
        [a, e, inc] = odst_elements(yf);
        ok = abs(a - r2) < 0.5 && e < 1e-4 && abs(inc - m.di_deg) < 1e-3 && abs(dvSim - o.dv) < 1e-4;
        report(ok, sprintf('%5.0f->%5.0f km di=%4.1f: %-38s a=%10.2f e=%.1e i=%7.3f dv_sim=%.5f dv=%.5f', ...
               m.alt1, m.alt2, m.di_deg, o.name, a, e, inc, dvSim, o.dv));
    end
end
fprintf('\n=== 3. Delta-v margin (low thrust needs 5.898 km/s for LEO -> GEO) ===\n');
mc = [6.0 0 1; 6.0 5 0; 6.5 5 1];
for i = 1:size(mc, 1)
    opts = odst_evaluate(odst_mission('dv_budget', mc(i, 1), 'margin_pct', mc(i, 2)));
    lt = opts(strcmp({opts.key}, 'low_thrust'));
    report(lt.feasible == logical(mc(i, 3)), sprintf('budget %.1f km/s, margin %g%%: low thrust feasible = %d', ...
           mc(i, 1), mc(i, 2), lt.feasible));
end
fprintf('\n=== 4. Gravity loss (Robbins), same values as python/validate.py ===\n');
r1 = RE + 400; at = (r1 + 42164)/2;
dv1 = sqrt(C.MU*(2/r1 - 1/at)) - sqrt(C.MU/r1);
check('LEO->GEO first burn loss, T/W 0.5 [m/s]', 1e3*odst_gravity_loss([r1 dv1], 320, 0.5), 14.94524, 1e-3);
check('LEO->GEO first burn loss, T/W 0.2 [m/s]', 1e3*odst_gravity_loss([r1 dv1], 320, 0.2), 93.40775, 1e-3);
opts = odst_evaluate(odst_mission('tw_chem', 0.2)); o = opts(strcmp({opts.key}, 'hohmann'));
check('Propellant fraction includes the gravity loss (Hohmann, T/W 0.2)', o.prop_frac, ...
      1 - exp(-(o.dv + o.dv_loss)/(o.isp*C.G0)), 1e-12);

fprintf('\n=== 4b. Edelbaum vs the averaged equations, integrated with Edelbaum''s yaw law ===\n');
f = 1e-6; r1 = RE + 400; r2 = 42164; di = 28.5*pi/180;         % 1e-6 km/s^2 = 1 mm/s^2
[~, tof] = odst_lowthrust(r1, r2, di, f*1e3);
V0 = sqrt(C.MU/r1); V1 = sqrt(C.MU/r2);
b0 = atan2(sin(pi/2*di), V0/V1 - cos(pi/2*di));
beta = @(t) atan2(V0*sin(b0), V0*cos(b0) - f*t);
[~, yy] = ode45(@(t, y) [-f*cos(beta(t)); 2/pi*f*sin(beta(t))/y(1)], [0 tof], [V0; 0], odeset('RelTol', 1e-11, 'AbsTol', 1e-13));
check('Edelbaum LEO -> GEO, di 28.5: speed after the closed-form TOF [km/s]', yy(end, 1), V1, 1e-8);
check('Edelbaum LEO -> GEO, di 28.5: plane change after the closed-form TOF [deg]', yy(end, 2)*180/pi, 28.5, 1e-6);
opts = odst_evaluate(odst_mission('di_deg', 120)); lt = opts(strcmp({opts.key}, 'low_thrust'));
report(~lt.applicable, sprintf('Edelbaum not offered beyond 114.6 deg (di 120): applicable = %d (%s)', lt.applicable, lt.note));

fprintf('\n=== 4c. Spacecraft model: masses, thrust, power, propellant on board ===\n');
% Thrust 0.5 m0 g0 is the generic T/W 0.5: every chemical option is unchanged, propellant = m0 x fraction
gen = odst_evaluate(odst_mission());
scm = odst_mission('mass0_kg', 1000, 'thrust_chem_N', 0.5*1000*C.G0*1e3);
sc = odst_evaluate(scm);
k = ~strcmp({gen.key}, 'low_thrust');
d = abs([gen(k).dv] - [sc(k).dv]) + abs([gen(k).dv_loss] - [sc(k).dv_loss]) + abs([gen(k).tof] - [sc(k).tof]);
check('Spacecraft with thrust = 0.5 m0 g0 vs generic T/W 0.5: chemical options [dv + loss + TOF difference]', max(d), 0, 1e-9);
check('  propellant = m0 x propellant fraction, Hohmann [kg]', sc(1).prop_kg, 1000*gen(1).prop_frac, 1e-9);
% Electric thrust from power, F = 2 eta P / ve; the spiral then runs at constant thrust as the mass falls.
% Integrate Edelbaum's averaged equations in time with a = F/m(t), yaw law in the velocity increment u.
F = odst_ep_thrust(scm);
check('Electric thrust 2 eta P / ve, 6.5 kW, eta 0.6, Isp 1600 s [N]', F, 2*0.6*6500/(1600*C.G0*1e3), 1e-12);
lt = sc(strcmp({sc.key}, 'low_thrust'));
r1 = RE + 400; r2 = RE + 35786; di = 28.5*pi/180; V0 = sqrt(C.MU/r1); V1 = sqrt(C.MU/r2);
b0 = atan2(sin(pi/2*di), V0/V1 - cos(pi/2*di)); ve = 1600*C.G0; f0 = F*1e-3/1000;   % km/s^2
yaw = @(u) atan2(V0*sin(b0), V0*cos(b0) - u);
rhs = @(~, y) [-f0/y(3)*cos(yaw(y(4))); 2/pi*f0/y(3)*sin(yaw(y(4)))/y(1); -f0/ve; f0/y(3)];
[~, yy] = ode45(rhs, [0 lt.tof], [V0; 0; 1; 0], odeset('RelTol', 1e-12, 'AbsTol', 1e-14));
check('Constant-thrust spiral LEO -> GEO, di 28.5: speed after the closed-form TOF [km/s]', yy(end, 1), V1, 1e-8);
check('  plane change after the closed-form TOF [deg]', yy(end, 2)*180/pi, 28.5, 1e-6);
check('  mass left, m/m0, vs exp(-dv/ve)', yy(end, 3), exp(-lt.dv/ve), 1e-9);
[~, tInf] = odst_lowthrust(r1, r2, di, 5e-4, 1e9); [~, tA] = odst_lowthrust(r1, r2, di, 5e-4);
check('Constant thrust with Isp -> infinity becomes constant acceleration [TOF ratio]', tInf/tA, 1, 1e-6);
% The propellant on board is a hard limit, at the same dv as the budget check (dv + loss, with the margin)
need = 1000*(1 - exp(-(sc(1).dv + sc(1).dv_loss)*1.05/(sc(1).isp*C.G0)));
on = odst_evaluate(odst_mission('mass0_kg', 1000, 'thrust_chem_N', 0.5*1000*C.G0*1e3, 'prop_chem_kg', need + 1e-6));
off = odst_evaluate(odst_mission('mass0_kg', 1000, 'thrust_chem_N', 0.5*1000*C.G0*1e3, 'prop_chem_kg', need - 1e-6));
report(on(1).feasible && ~off(1).feasible, sprintf(['Hohmann needs %.6f kg with the 5%% margin: feasible with ' ...
       '1e-6 kg more on board, not with 1e-6 kg less'], need));
[t, rec] = odst_evaluate(odst_mission('mass0_kg', 1000, 'prop_chem_kg', 600, 'prop_ep_kg', 350));
report(strcmp(rec.key, 'low_thrust') && ~any([t(~strcmp({t.key}, 'low_thrust')).feasible]), sprintf(['600 kg of ' ...
       'chemical propellant on board: every chemical option is out (%s); low thrust is recommended, %.0f kg in %.1f days'], ...
       t(1).reasons, rec.prop_kg, rec.tof/86400));
% Same as python/odst.py, default spacecraft, LEO -> GEO: [dv km/s, loss km/s, TOF s, propellant kg]
pyRefSc = [4.221535383 0.014972551 19048.5625 740.762642; 4.19629709 0.015327503 19048.5625 738.698885
           5.470508212 0.01867395 10862.2498 826.085245; 5.897522253 0 9888981.43 313.303212];
o = odst_evaluate(odst_mission('mass0_kg', 1000)); o = o([o.applicable]);
got = [[o.dv].', [o.dv_loss].', [o.tof].', [o.prop_kg].'];
report(isequal(size(got), size(pyRefSc)) && max(max(abs(got(:, [1 2 4]) - pyRefSc(:, [1 2 4])))) < 1e-6 && ...
       max(abs(got(:, 3) - pyRefSc(:, 3))) < 1e-3, sprintf(['Default spacecraft, LEO -> GEO: dv, loss, TOF and ' ...
       'propellant of %d options as Python (largest propellant difference %.1e kg)'], size(got, 1), ...
       max(abs(got(:, 4) - pyRefSc(:, 4)))));
% The export carries the spacecraft: GMAT draws each burn's propellant from a tank (python/gmat_crosscheck.py)
txt = odst_gmat_script(scm, sc(1)); js = jsondecode(odst_burns_json(scm, sc(1)));
report(contains(txt, 'Sat.DryMass = 200;') && contains(txt, 'Fuel.FuelMass = 800;') && ...
       count(txt, 'DecrementMass = true;') == 2 && js.spacecraft.mass0_kg == 1000 && js.spacecraft.propellant_kg == sc(1).prop_kg, ...
       sprintf(['GMAT script with the spacecraft''s dry mass, an 800 kg tank and mass-decrementing burns; JSON with ' ...
                'the propellant, %.1f kg'], sc(1).prop_kg));

fprintf('\n=== 4d. Uncertainty: Monte Carlo over Isp, thrust and pointing errors; the Pareto front ===\n');
u = odst_mrg32k3a(5, [0 0 1234567 0 0 1234567]);
report(isequal(round(u.'/2.328306549295727688e-10) + 1, [1459213977 2827710106 4245671317 3877608661 2595287583]), ...
       'MRG32k3a reproduces the published test vector (seed 1234567, first five outputs)');
% the same normals as Python: first and last run
Z = odst_mc_normals();
py = [-0.847924823347079 1.8460727873862615 0.7028567229701457 -1.361475967116543 -1.6978660974898043 -0.3615845453046289 -0.7451327634504835 0.953673039659349
      0.3355242415679864 0.630794942430949 -1.095285273790694 1.0227744211709222 0.7098950824471811 0.858224862332908 -0.7171417966521667 -0.10778190475429406];
check('4000 runs: the first 8 normals of the first and last run vs Python [largest difference]', max(max(abs(Z([1 end], 1:8) - py))), 0, 1e-12);
check('  sum of the first 8 normals of every run (32000) vs Python', sum(sum(Z(:, 1:8))), -224.82505288590977, 1e-9);
% with every error zero, each run is the nominal one
zero = {'uncertainty', true, 'sigma_isp_pct', 0, 'sigma_thrust_pct', 0, 'sigma_point_deg', 0};
ms = {{}, {'mass0_kg', 1000, 'prop_chem_kg', 500, 'prop_ep_kg', 350}, {'alt1', 250, 'apo1', 35786, 'propulsion', 'chemical'}, ...
      {'weights', [0.25 0.55 0.05 0.15], 'tof_max_h', 3300}};
ok = true;
for c = 1:numel(ms)
    a = odst_evaluate(odst_mission(ms{c}{:})); b = odst_evaluate(odst_mission(ms{c}{:}, zero{:}));
    ok = ok && isequal([a.feasible], [b.feasible]) && all([b([b.feasible]).chance] == 1);
end
report(ok, 'Errors set to zero, 4 missions: every feasible option fits in all runs, and feasibility is unchanged');
% the same chances as Python, run for run (runs that fit, of 4000)
mcCases = {'Spacecraft near its tank limits', {'mass0_kg', 1000, 'prop_chem_kg', 765, 'prop_ep_kg', 330, 'tof_max_h', 2800, 'confidence_pct', 0}, ...
              {'hohmann', 3328; 'hohmann_split', 3740; 'low_thrust', 3078}
           'Propellant-critical, 3300 h', {'weights', [0.25 0.55 0.05 0.15], 'tof_max_h', 3300}, ...
              {'hohmann', 4000; 'hohmann_split', 4000; 'one_tangent', 4000; 'low_thrust', 2577}
           'Small tank, 2800 h', {'mass0_kg', 1000, 'prop_chem_kg', 500, 'prop_ep_kg', 350, 'tof_max_h', 2800}, {'low_thrust', 3325}
           'GTO -> GEO, budget 1.95 km/s', {'alt1', 250, 'apo1', 35786, 'propulsion', 'chemical', 'dv_budget', 1.95}, {'apse', 2725; 'lambert', 2708}
           'LEO -> GEO, budget 4.5 km/s, pointing 1.5 deg', {'dv_budget', 4.5, 'sigma_point_deg', 1.5}, {'hohmann', 114; 'hohmann_split', 440}};
for c = 1:size(mcCases, 1)
    opts = odst_evaluate(odst_mission('uncertainty', true, mcCases{c, 2}{:}));
    got = opts(~isnan([opts.chance])); ref = mcCases{c, 3};
    txt = strjoin(arrayfun(@(o) sprintf('%s %d', o.key, round(o.chance*4000)), got, 'UniformOutput', false), ', ');
    report(isequal({got.key}, ref(:, 1).') && isequal(round([got.chance]*4000), [ref{:, 2}]), ...
           sprintf('%s: runs that fit %s (Python %s)', mcCases{c, 1}, txt, ...
                   strjoin(cellfun(@(k, n) sprintf('%s %d', k, n), ref(:, 1).', ref(:, 2).', 'UniformOutput', false), ', ')));
end
% the rule: the spiral makes the 3300 h limit in 64 % of runs, so it is out and Hohmann is recommended
[~, a] = odst_evaluate(odst_mission('weights', [0.25 0.55 0.05 0.15], 'tof_max_h', 3300));
[opts, b] = odst_evaluate(odst_mission('weights', [0.25 0.55 0.05 0.15], 'tof_max_h', 3300, 'uncertainty', true));
lt = opts(strcmp({opts.key}, 'low_thrust'));
report(strcmp(a.key, 'low_thrust') && ~lt.feasible && strcmp(b.key, 'hohmann') && ...
       strcmp(lt.reasons, 'fits the limits in 64% of runs with the errors, 90% required'), ...
       sprintf('Propellant-critical, 3300 h: %s without errors; with them: %s; recommended: %s', a.name, lt.reasons, b.name));
opts = odst_evaluate(odst_mission());
front = {opts([opts.pareto]).key};
report(isequal(front, {'hohmann_split', 'one_tangent', 'low_thrust'}), ['Pareto front of LEO -> GEO: ' strjoin(front, ', ')]);

fprintf('\n=== 4e. J2: the planes'' drift, the spiral, the nodal drift, against Python''s values ===\n');
PC = [0.25 0.55 0.05 0.15];
% name, mission, recommended, {key, dv, tof, applicable, J2 dv, parking altitude, wait days, corrections} from odst.py
j2Cases = { ...
    'raan10', {'alt1', 500, 'alt2', 500, 'elements', [51.6 0 0 51.6 10 0], 'tof_max_h', 384, 'weights', PC}, 'nodal_drift', {{'hohmann', 1.0399338192964325, 0.0, true, 0.0, NaN, NaN, NaN}; ...
          {'hohmann_split', 1.0399338192964325, 0.0, false, NaN, NaN, NaN, NaN}; ...
          {'bielliptic', 1.0399338223227859, 5676.9780327835915, false, NaN, NaN, NaN, NaN}; ...
          {'low_thrust', 2.4318678440890706, 4863735.688178142, true, 0.8002103036610191, NaN, NaN, 2.0}; ...
          {'nodal_drift', 0.3042482984417413, 1382401.3823999988, true, NaN, 783.388963230549, 15.9322694368591, NaN}}; ...
    'raan10inc', {'alt1', 500, 'alt2', 500, 'elements', [51.6 0 0 53 10 0], 'tof_max_h', 384, 'weights', PC}, 'nodal_drift', {{'hohmann', 1.066152021082588, 0.0, true, 0.0, NaN, NaN, NaN}; ...
          {'hohmann_split', 1.066152021082588, 0.0, false, NaN, NaN, NaN, NaN}; ...
          {'bielliptic', 1.0661520240892794, 5676.9780327835915, false, NaN, NaN, NaN, NaN}; ...
          {'low_thrust', 1.9929670917756679, 3985934.183551336, true, 0.32027102790091955, NaN, NaN, 2.0}; ...
          {'nodal_drift', 0.49060539494954797, 1382401.3823999998, true, NaN, 858.484512771689, 15.931725164999492, NaN}}; ...
    'spiral', {'alt1', 400, 'alt2', 1200, 'elements', [51.6 0 0 51.6 0 0], 'tof_max_h', 1e5}, 'hohmann', {{'hohmann', 0.4157515124549431, 3026.206774746056, true, 1.5377694229656846e-05, NaN, NaN, NaN}; ...
          {'hohmann_split', 0.4157361347607134, 3026.206774746056, false, NaN, NaN, NaN, NaN}; ...
          {'bielliptic', 0.4157361383783885, 6308.857407789015, false, NaN, NaN, NaN, NaN}; ...
          {'one_tangent', 2.4961364736241793, 916.4530230408586, true, 6.531710199375595e-07, NaN, NaN, NaN}; ...
          {'low_thrust', 1.7514522899225888, 3502904.579845178, true, 1.335392849131551, NaN, NaN, 3.0}; ...
          {'nodal_drift', NaN, NaN, false, NaN, NaN, NaN, NaN}}; ...
    'spiralsc', {'alt1', 400, 'alt2', 1200, 'elements', [51.6 0 0 51.6 0 0], 'tof_max_h', 1e5, 'mass0_kg', 1000, 'prop_ep_kg', 150}, 'hohmann', {{'hohmann', 0.4157515124549431, 3026.206774746056, true, 1.5377694229656846e-05, NaN, NaN, NaN}; ...
          {'hohmann_split', 0.4157361347607134, 3026.206774746056, false, NaN, NaN, NaN, NaN}; ...
          {'bielliptic', 0.4157361383783885, 6308.857407789015, false, NaN, NaN, NaN, NaN}; ...
          {'one_tangent', 2.4961364736241793, 916.4530230408586, true, 6.531710199375595e-07, NaN, NaN, NaN}; ...
          {'low_thrust', 1.736063744045925, 3306034.2183484803, true, 1.3200043032548872, NaN, NaN, 3.0}; ...
          {'nodal_drift', NaN, NaN, false, NaN, NaN, NaN, NaN}}; ...
    'slow', {'alt1', 400, 'alt2', 1200, 'elements', [51.6 0 0 51.6 0 0], 'tof_max_h', 1e5, 'accel_ep', 1e-4}, 'hohmann', {{'hohmann', 0.4157515124549431, 3026.206774746056, true, 1.5377694229656846e-05, NaN, NaN, NaN}; ...
          {'hohmann_split', 0.4157361347607134, 3026.206774746056, false, NaN, NaN, NaN, NaN}; ...
          {'bielliptic', 0.4157361383783885, 6308.857407789015, false, NaN, NaN, NaN, NaN}; ...
          {'one_tangent', 2.4961364736241793, 916.4530230408586, true, 6.531710199375595e-07, NaN, NaN, NaN}; ...
          {'low_thrust', NaN, NaN, false, NaN, NaN, NaN, NaN}; ...
          {'nodal_drift', NaN, NaN, false, NaN, NaN, NaN, NaN}}; ...
    'meo', {'alt1', 400, 'alt2', 20200, 'elements', [28.5 0 0 55 40 0], 'tof_max_h', 1e5}, 'hohmann', {{'hohmann', 4.406918154306862, 10717.717121050175, true, 0.001485716537419357, NaN, NaN, NaN}; ...
          {'hohmann_split', 4.330759229284638, 10717.717121050175, true, 0.0014160622915824206, NaN, NaN, NaN}; ...
          {'bielliptic', 4.405432438024526, 32278.66198453821, false, NaN, NaN, NaN, NaN}; ...
          {'one_tangent', 5.820344649743336, 6055.36761364423, true, 0.00029958441251221046, NaN, NaN, NaN}; ...
          {'low_thrust', 14.719383419439804, 29438766.83887961, true, 8.254497599477219, NaN, NaN, 3.0}; ...
          {'nodal_drift', 4.006970399058953, 3937656.8604975925, true, NaN, 400.0, 45.45068452982109, NaN}}; ...
    'descent', {'alt1', 20200, 'alt2', 500, 'elements', [55 0 0 55 3 0], 'tof_max_h', 1e5}, 'hohmann', {{'hohmann', 3.376835683771683, 10765.949782789343, true, -0.0019006893363155974, NaN, NaN, NaN}; ...
          {'hohmann_split', 3.3761343495235128, 10765.949782789343, true, -0.0016026078272135535, NaN, NaN, NaN}; ...
          {'bielliptic', 3.3787363743786605, 32326.89464629656, false, NaN, NaN, NaN, NaN}; ...
          {'one_tangent', 6.742805504326066, 8981.535790465627, true, -0.0009324085572508523, NaN, NaN, NaN}; ...
          {'low_thrust', NaN, NaN, false, NaN, NaN, NaN, NaN}; ...
          {'nodal_drift', 3.3723859323697707, 60628.84854822148, true, NaN, 20200.0, 0.5771168838591683, NaN}}; ...
    'geo', {'alt1', 400, 'alt2', 35786, 'elements', [51.6 120 0 0 0 0]}, 'hohmann', {{'hohmann', 4.82464984429949, 19048.562509797113, true, 0.0, NaN, NaN, NaN}; ...
          {'hohmann_split', 4.784537299351257, 19048.562509797113, true, 0.0, NaN, NaN, NaN}; ...
          {'bielliptic', 4.501603227731081, 920556.7324006476, true, 0.0, NaN, NaN, NaN}; ...
          {'one_tangent', 5.901213052479644, 10862.249834254219, true, 0.0, NaN, NaN, NaN}; ...
          {'low_thrust', 7.805548981621778, 15611097.963243555, true, 0.0, NaN, NaN, 0.0}; ...
          {'nodal_drift', NaN, NaN, false, NaN, NaN, NaN, NaN}}; ...
    'kourou', {'alt1', 250, 'apo1', 35786, 'alt2', 35786, 'elements', [6 0 178 0 0 0], 'propulsion', 'chemical'}, 'lambert', {{'apse', NaN, NaN, false, NaN, NaN, NaN, NaN}; ...
          {'lambert', 1.4967607162136738, 3479.8680687829756, true, 0.0, NaN, NaN, NaN}}};
close = @(a, b, tol) (isnan(a) && isnan(b)) || abs(a - b) <= tol*max(1, abs(b));
for c = 1:size(j2Cases, 1)
    [opts, rec] = odst_evaluate(odst_mission('j2', true, j2Cases{c, 2}{:}));
    want = j2Cases{c, 4}; bad = {};
    for q = 1:size(want, 1)
        w = want{q}; o = opts(strcmp({opts.key}, w{1}));
        if isempty(o), bad{end+1} = [w{1} ' missing']; continue; end %#ok<AGROW>
        if o.applicable ~= w{4}, bad{end+1} = sprintf('%s applicable %d', w{1}, o.applicable); end %#ok<AGROW>
        if ~w{4}, continue; end
        if ~close(o.dv, w{2}, 1e-9) || ~close(o.tof, w{3}, 1e-9), bad{end+1} = sprintf('%s dv %.12g tof %.10g', w{1}, o.dv, o.tof); end %#ok<AGROW>
        if ~isnan(w{5}) && ~close(o.j2_dv, w{5}, 1e-9), bad{end+1} = sprintf('%s J2 dv %.12g', w{1}, o.j2_dv); end %#ok<AGROW>
        if ~isnan(w{6}) && (~close(o.params.alt_d, w{6}, 1e-9) || ~close(o.params.wait_days, w{7}, 1e-9))
            bad{end+1} = sprintf('%s parking %.9g km %.9g d', w{1}, o.params.alt_d, o.params.wait_days); %#ok<AGROW>
        end
        if ~isnan(w{8}) && o.params.corrections ~= w{8}, bad{end+1} = sprintf('%s corrections %d', w{1}, o.params.corrections); end %#ok<AGROW>
    end
    rk = ''; if ~isempty(rec), rk = rec.key; end
    if ~strcmp(rk, j2Cases{c, 3}), bad{end+1} = ['recommended ' rk]; end %#ok<AGROW>
    report(isempty(bad), sprintf('J2, %s: %d options as Python (%s recommended)%s', j2Cases{c, 1}, numel(opts), rk, ...
           strjoin(strcat({' '}, bad), ';')));
end
% from a polar orbit the parking orbits tie on time: the cheapest is taken (Python: parking where it starts)
opts = odst_evaluate(odst_mission('j2', true, 'alt1', 452, 'alt2', 1273, 'elements', [90 6 0 84.2 235.2 0], 'tof_max_h', 20));
o = opts(strcmp({opts.key}, 'nodal_drift'));
report(o.applicable && ~o.feasible && close(o.dv, 0.9626656591468357, 1e-9) && close(o.tof, 37182873.32233598, 1e-9) && close(o.params.alt_d, 452, 1e-9), ...
       sprintf('J2, polar start, time limit out of reach: nodal drift %.4f km/s parking at %.0f km, as Python', o.dv, o.params.alt_d));
% the local burns J2 plans fly in: VNB components, events (0 time, 1 apoapsis, 2 periapsis, 3 node), the node-signed burn
events = {'time', 'apoapsis', 'periapsis', 'node'};
lbCases = {'raan10inc', 'nodal_drift', [0.0, 0.09606579681241367, 0.0, 0.0, 0, 0; 2950.1640720218165, 0.09485322144623076, 0.0, 0.0, 1, 0; 1379451.218327978, -0.0970403476500783, 0.17900980832924193, 0.0, 3, 1; 1382401.3823999998, -0.09606579681241278, 0.0, 0.0, 2, 0]; 'descent', 'hohmann', [0.0, -1.3909970758888304, 0.08910037861377262, -6.938893903907228e-18, 0, 0; 10765.949782789343, -1.982987864746827, -2.7755575615628914e-17, 6.938893903907228e-18, 2, 0]};
for c = 1:2
    args = j2Cases{strcmp(j2Cases(:, 1), lbCases{c, 1}), 2};
    m = odst_mission('j2', true, args{:}); opts = odst_evaluate(m); o = opts(strcmp({opts.key}, lbCases{c, 2}));
    [~, b] = odst_local_burns(m, o); w = lbCases{c, 3};
    ok = numel(b) == size(w, 1);
    for q = 1:min(numel(b), size(w, 1))
        ok = ok && abs(b(q).t - w(q, 1)) <= 1e-6*max(1, w(q, 1)) && max(abs(b(q).dv - w(q, 2:4))) < 1e-12 && ...
             strcmp(b(q).at, events{w(q, 5) + 1}) && b(q).signed == logical(w(q, 6));
    end
    report(ok, sprintf('J2, %s %s: local burns as Python (%s)', lbCases{c, 1}, lbCases{c, 2}, strjoin({b.at}, ', ')));
end
% the export: a J2 field, VNB burns, the node crossing, the target alongside
args = j2Cases{2, 2}; m = odst_mission('j2', true, args{:}); opts = odst_evaluate(m); o = opts(strcmp({opts.key}, 'nodal_drift'));
txt = odst_gmat_script(m, o); js = jsondecode(odst_burns_json(m, o));
bb = js.burns; if ~iscell(bb), bb = num2cell(bb); end      % a cell array when the burns' fields differ
report(contains(txt, 'GravityField.Earth.Degree = 2;') && count(txt, 'Axes = VNB;') == 5 && contains(txt, '{Sat.EarthMJ2000Eq.Z = 0}') && ...
       contains(txt, 'If Sat.EarthMJ2000Eq.VZ > 0') && contains(txt, 'Create Spacecraft TargetSat;') && ...
       isequal(cellfun(@(x) x.at, bb(:).', 'UniformOutput', false), {'time', 'apoapsis', 'node', 'periapsis'}), ...
       'J2 export of the nodal drift: J2 field, 5 VNB burns, the node crossing, TargetSat; JSON events time, apoapsis, node, periapsis');
% J2 needs absolute inclinations: with only an inclination change it changes nothing
a = odst_evaluate(odst_mission('alt2', 1200, 'di_deg', 4, 'j2', true)); b = odst_evaluate(odst_mission('alt2', 1200, 'di_deg', 4));
report(isequal({a.key}, {b.key}) && isequal([a.dv], [b.dv]), 'J2 on with only an inclination change: the two-body options, unchanged');

fprintf('\n=== 5. Elliptical orbits and Lambert ===\n');
% Curtis, Orbital Mechanics for Engineering Students, Example 5.2
r1v = [5000; 10000; 2100]; r2v = [-14600; 2500; 7000];
[v1, v2] = odst_lambert(r1v, r2v, 3600);
ref = [-5.9925 1.9254 3.2456 -3.3125 -4.1966 -0.38529]; got = [v1; v2];
for k = 1:6
    check(sprintf('Lambert, Curtis Ex. 5.2, v%d(%d) [km/s]', 1 + (k > 3), mod(k - 1, 3) + 1), got(k), ref(k), 1e-4);
end
[~, yy] = ode45(@(~, y) [y(4:6); -C.MU*y(1:3)/norm(y(1:3))^3], [0 3600], [r1v; v1], ...
                odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
check('Lambert v1 propagated for 1 h reaches r2 [km]', norm(yy(end, 1:3).' - r2v), 0, 1e-6);
di = 20*pi/180; ap = odst_apse(7000, 7000, 14000, 14000, di, 0);
check('Apse-to-apse between circular orbits == Hohmann, optimal split [km/s]', ap.dv, ...
      odst_hohmann(7000, 14000, di, odst_hohmann_split_opt(7000, 14000, di)), 1e-9);
lo = odst_lambert_opt(7000, 7000, 14000, 14000, 0, 0, 1e9);
check('Lambert optimum, coplanar 7000 -> 14000 km, vs Hohmann [m/s]', 1e3*lo.dv, 1e3*odst_hohmann(7000, 14000, 0), 1.0);
% Same missions as python/validate.py; reference dv [km/s] and TOF [s] from python/odst.py
ell = {odst_mission('alt1', 250, 'apo1', 35786, 'alt2', 35786, 'di_deg', 28.5), ...
       odst_mission('alt1', 400, 'alt2', 600, 'apo2', 39700, 'dw_deg', 270, 'di_deg', 0), ...
       odst_mission('alt1', 500, 'apo1', 5000, 'alt2', 800, 'apo2', 12000, 'dw_deg', 60, 'di_deg', 0), ...
       odst_mission('alt1', 500, 'apo1', 5000, 'alt2', 800, 'apo2', 12000, 'dw_deg', 60, 'di_deg', 10), ...
       odst_mission('alt1', 300, 'apo1', 12500, 'alt2', 350, 'apo2', 8600, 'dw_deg', 90, 'di_deg', 0, ...
                    'tof_max_h', 1950/3600)};     % 100 km floor active
pyRef = {[1.833344084 0; 1.833344097 20504.3564], [2.476172306 21378.6762; 2.476173769 21646.7968], ...
         [NaN NaN; 1.183670980 9359.7970], [NaN NaN; 1.859013203 3601.2428], [NaN NaN; 2.546787403 1629.7306]};
pyRec = {'apse', 'apse', 'lambert', 'lambert', 'lambert'};
tolT = [1 1 1 1 5];      % s; on the 100 km floor the optimum's time is less sharply determined
for c = 1:numel(ell)
    m = ell{c};
    [opts, rec] = odst_evaluate(m);
    mt = sprintf('%.0fx%.0f -> %.0fx%.0f km di=%g dw=%g', m.alt1, max(m.alt1, m.apo1), ...
                  m.alt2, max(m.alt2, m.apo2), m.di_deg, m.dw_deg);
    report(strcmp(rec.key, pyRec{c}), sprintf('%s: recommends %s (Python: %s)', mt, rec.key, pyRec{c}));
    for i = 1:numel(opts)
        o = opts(i);
        if ~o.applicable
            report(isnan(pyRef{c}(i, 1)), sprintf('%s: %s not applicable (%s)', mt, o.key, o.note));
            continue
        end
        report(abs(o.dv - pyRef{c}(i, 1)) < 1e-6 && abs(o.tof - pyRef{c}(i, 2)) < tolT(c), ...
               sprintf('%s: %-8s dv %.9f vs Python %.9f (%.1e km/s), TOF %.4f vs %.4f s', mt, o.key, ...
                       o.dv, pyRef{c}(i, 1), abs(o.dv - pyRef{c}(i, 1)), o.tof, pyRef{c}(i, 2)));
        [~, ~, yf, dvSim] = odst_simulate_vec(o.y0, o.vplan, 3600);
        [dh, de] = targetMismatch(yf, m);
        report(dh < 1e-6 && de < 1e-6 && abs(dvSim - o.dv) < 1e-9, ...
               sprintf('%s: %-8s propagated arrival |dh|/h = %.1e, |de| = %.1e', mt, o.key, dh, de));
        if strcmp(o.key, 'lambert')            % the arc must stay above the 100 km floor
            [~, ~, y2] = odst_simulate_vec([o.y0(1:3); o.y0(4:6) + o.vplan(1, 2:4).'], zeros(0, 4), o.vplan(2, 1));
            low = odst_arc_min_radius(o.y0(1:3), o.y0(4:6) + o.vplan(1, 2:4).', y2(1:3)) - RE;
            report(low >= 100 - 1e-6, sprintf('%s: lambert   arc lowest altitude %.4f km (floor 100 km)', mt, low));
        end
    end
end
fprintf('\n=== 5b. Export: exact burn plans, GMAT script and JSON ===\n');
% Kepler propagation vs Kepler's equation, start and end points from eccentric anomalies [E0 E1]
arcs = [6778 6778 0.3 2.0; 7000 42164 1.0 4.5; 6678 384400 5.0 7.0; 8000 30000 2.5 3.0];   % rp ra [km]
worst = 0;
for q = 1:size(arcs, 1)
    rp = arcs(q, 1); ra = arcs(q, 2); a = (rp + ra)/2; e = (ra - rp)/(ra + rp);
    th = @(E) 2*atan2(sqrt(1 + e)*sin(E/2), sqrt(1 - e)*cos(E/2));
    t = @(E) sqrt(a^3/C.MU)*(E - e*sin(E));
    [r0, v0] = odst_orbit_state(rp, ra, 0.4, 1.1, th(arcs(q, 3)));
    y = odst_kepler([r0; v0], t(arcs(q, 4)) - t(arcs(q, 3)));
    [rv, vv] = odst_orbit_state(rp, ra, 0.4, 1.1, th(arcs(q, 4)));
    worst = max([worst, norm(y(1:3) - rv)/norm(rv), norm(y(4:6) - vv)/norm(vv)]);
end
check('Kepler propagation vs Kepler''s equation, 4 arcs up to e = 0.97 [1e-12 relative]', 1e12*worst, 0, 100);
% ... and hyperbolic arcs (fast Lambert transfers), from the hyperbolic anomaly: M = e sinh F - F; [rp e F0 F1]
arcs = [7000 1.05 -1 2; 12000 3 -2.5 0.5; 20000 7.5 0.2 2.9];
worst = 0;
for q = 1:size(arcs, 1)
    rp = arcs(q, 1); e = arcs(q, 2); a = rp/(1 - e); p = rp*(1 + e);
    [P, Q] = odst_perifocal(0.4, 1.1, 0.7);
    th = @(F) 2*atan(sqrt((e + 1)/(e - 1))*tanh(F/2));
    t = @(F) sqrt((-a)^3/C.MU)*(e*sinh(F) - F);
    st = @(nu) [p/(1 + e*cos(nu))*(cos(nu)*P + sin(nu)*Q); sqrt(C.MU/p)*(-sin(nu)*P + (e + cos(nu))*Q)];
    y = odst_kepler(st(th(arcs(q, 3))), t(arcs(q, 4)) - t(arcs(q, 3))); ref = st(th(arcs(q, 4)));
    worst = max([worst, norm(y(1:3) - ref(1:3))/norm(ref(1:3)), norm(y(4:6) - ref(4:6))/norm(ref(4:6))]);
end
check('Kepler propagation vs the hyperbolic Kepler equation, 3 arcs, e 1.05 to 7.5 [1e-12 relative]', 1e12*worst, 0, 100);
% The circular models' exact vector plans: propagated, they reach the target plane itself (normal
% [0; -sin(di); cos(di)]), not only its inclination, at the closed-form cost
more = {odst_mission('alt1', 400, 'alt2', 400, 'di_deg', 60), ...
        odst_mission('alt1', 35786, 'alt2', 400, 'di_deg', 28.5, 'tof_max_h', 3)};
allc = [cases, more];
for c = 1:numel(allc)
    m = allc{c};
    opts = odst_evaluate(m);
    r2 = RE + m.alt2;
    for i = 1:numel(opts)
        o = opts(i);
        if isempty(o.plan) || ~o.applicable, continue; end
        [~, ~, yf, dvSim] = odst_simulate_vec(o.y0, o.vplan, 3600);
        [dh, de] = targetMismatch(yf, m);
        report(dh < 1e-6 && de < 1e-6 && abs(dvSim - o.dv) < 1e-12, ...
               sprintf('%5.0f->%5.0f km di=%4.1f: %-38s vector plan |dh|/h=%.1e |de|=%.1e, sum|dv| - dv=%.1e', ...
                       m.alt1, m.alt2, m.di_deg, o.name, dh, de, dvSim - o.dv));
    end
end
% Same plans as python/odst.py: closed-form burns (no optimiser), start state and rows [t dvx dvy dvz]
opts = odst_evaluate(cases{1});
pyRefV = {'hohmann', [6778.137 0 0 0 7.668558175407055 0], ...
          [0 0 2.3974725168152418 0; 19048.562509797113 9.211582468598876e-16 -1.0838903670048219 -1.467101568882389]
          'one_tangent', [-5964.52423775776 -3219.874501897429 0 3.6428586402591523 -6.748063827253831 0], ...
          [0 1.2532385966758652 -2.3215103511524484 0; 10862.249834254219 -2.18568965499367 1.0451944822253934 1.467101568882389]};
for q = 1:size(pyRefV, 1)
    o = opts(strcmp({opts.key}, pyRefV{q, 1}));
    d = max([abs(o.y0(:).' - pyRefV{q, 2})./max(1, abs(pyRefV{q, 2})), ...
             reshape(abs(o.vplan - pyRefV{q, 3})./max(1, abs(pyRefV{q, 3})), 1, [])]);
    report(d < 1e-12, sprintf('LEO -> GEO %-12s start state and vector plan vs Python: largest difference %.1e', o.key, d));
end
% GMAT script and JSON read back to the exact plan (GMAT flies them: python/gmat_crosscheck.py)
for mm = {cases{1}, ell{2}, odst_mission('alt1', 400, 'alt2', 600, 'apo2', 39700, 'elements', [51.6 30 0 63.4 100 270], 'tof_max_h', 30)}
    m = mm{1};
    opts = odst_evaluate(m);
    for i = 1:numel(opts)
        o = opts(i);
        if isempty(o.vplan), continue; end
        txt = odst_gmat_script(m, o); js = jsondecode(odst_burns_json(m, o));
        [y0i, plani] = odst_inertial_plan(m, o);
        val = @(k) str2double(regexp(txt, ['^' regexptranslate('escape', k) ' = (\S+);$'], 'tokens', 'once', 'lineanchors'));
        lines = splitlines(txt(strfind(txt, 'BeginMissionSequence'):end));
        P = zeros(0, 4); t = 0;
        for q = 1:numel(lines)
            if startsWith(lines{q}, 'Propagate') && size(P, 1) < size(plani, 1)
                t = t + str2double(regexp(lines{q}, 'ElapsedSecs = (\S+)}', 'tokens', 'once'));
            elseif startsWith(lines{q}, 'Maneuver')
                b = size(P, 1) + 1;
                P(b, :) = [t, val(sprintf('Burn%d.Element1', b)), val(sprintf('Burn%d.Element2', b)), ...
                           val(sprintf('Burn%d.Element3', b))];
            end
        end
        y0 = cellfun(@(k) val(['Sat.' k]), {'X', 'Y', 'Z', 'VX', 'VY', 'VZ'}).';
        jb = js.burns;
        ok = val('Earth.Mu') == C.MU && isequal(y0, y0i) && isequal(js.start_state_km_kms, y0i) && ...
             isequal(size(P), size(plani)) && all(abs(P(:, 1) - plani(:, 1)) <= 1e-9*max(1, plani(:, 1))) && ...
             isequal(P(:, 2:4), plani(:, 2:4)) && isequal([[jb.t_s].', [jb.dv_kms].'], plani) && ...
             (isempty(m.elements) && isempty(js.elements_deg) || isequal(js.elements_deg(:).', m.elements));
        report(ok, sprintf('%-38s GMAT script and JSON give back the plan (%d burns, start state, Earth mu%s)', ...
                           o.name, size(P, 1), repmat(', elements', 1, ~isempty(m.elements))));
    end
end

fprintf('\n=== 5c. Orbital elements (inclination, RAAN, argument of periapsis) ===\n');
dvs = @(m) evalDv(odst_evaluate(m));
% Elements that put the start orbit in x-y with the node on +x are the relative inputs themselves
rel = {odst_mission('alt1', 400, 'alt2', 600, 'apo2', 39700, 'di_deg', 20, 'dw_deg', 270), [0 0 0 20 0 270]
       odst_mission('alt1', 250, 'apo1', 35786, 'alt2', 35786, 'di_deg', 28.5), [0 0 0 28.5 0 0]};
for q = 1:size(rel, 1)
    r = rel{q, 1};
    m = odst_mission('alt1', r.alt1, 'apo1', r.apo1, 'alt2', r.alt2, 'apo2', r.apo2, 'elements', rel{q, 2});
    check(sprintf('Elements [%s] == relative inputs, every option [km/s]', strtrim(sprintf('%g ', rel{q, 2}))), ...
          max(abs(dvs(m) - dvs(r))), 0, 1e-12);
end
% Turning both orbits together (both RAANs + 40 deg) changes nothing
c = dvs(odst_mission('alt1', 300, 'apo1', 12000, 'alt2', 800, 'apo2', 20000, 'elements', [28.5 10 30 40 70 120], 'tof_max_h', 20));
d = dvs(odst_mission('alt1', 300, 'apo1', 12000, 'alt2', 800, 'apo2', 20000, 'elements', [28.5 50 30 40 110 120], 'tof_max_h', 20));
check('Both orbits turned 40 deg about the pole: Lambert, 3D [km/s]', max(abs(c - d)), 0, 1e-9);
% A RAAN change alone: the planes are 2 asin(sin i sin(dRAAN/2)) apart (spherical trigonometry)
m = odst_mission('alt1', 500, 'alt2', 500, 'elements', [51.6 0 0 51.6 10 0]);
ref = 2*asin(sin(51.6*pi/180)*sin(5*pi/180));
check('RAAN 0 -> 10 deg at i = 51.6 deg: angle between the planes [deg]', m.di_deg, ref*180/pi, 1e-12);
opts = odst_evaluate(m);
check('  its single-burn plane change, 2 v sin(di/2) [km/s]', opts(1).dv, 2*sqrt(C.MU/(RE + 500))*sin(ref/2), 1e-12);
% GTO with its apogee on the node line equals the relative GTO; 2 deg off it, only Lambert is exact
g = odst_evaluate(odst_mission('alt1', 250, 'apo1', 35786, 'alt2', 35786, 'elements', [6 0 180 0 0 0]));
g6 = odst_evaluate(odst_mission('alt1', 250, 'apo1', 35786, 'alt2', 35786, 'di_deg', 6));
check('GTO i 6, argp 180 -> GEO: apse-to-apse == relative di 6 [km/s]', g(1).dv, g6(1).dv, 1e-12);
o = odst_evaluate(odst_mission('alt1', 250, 'apo1', 35786, 'alt2', 35786, 'elements', [6 0 178 0 0 0]));
report(~o(1).applicable && o(2).applicable && o(2).dv > g(1).dv, sprintf(['GTO i 6, argp 178 -> GEO: apse-to-apse ' ...
       'not offered (%s); Lambert %.6f km/s, %.1f m/s above the aligned case'], o(1).note, o(2).dv, 1e3*(o(2).dv - g(1).dv)));
% Same derived geometry and results as python/odst.py: [di w1 dw] deg, then [dv km/s, TOF s] per applicable option
pyEl = {{'alt1', 400, 'alt2', 35786, 'elements', [51.6 120 0 0 0 0]}, [51.60000000000001 0 0], ...
        [4.824649844 19048.5625; 4.784537299 19048.5625; 4.501603228 920556.7324; 5.901213052 10862.2498; 7.805548982 15611097.9632]
        {'alt1', 250, 'apo1', 35786, 'alt2', 35786, 'elements', [6 0 178 0 0 0]}, [6.000000000000001 358 0], [1.496760716 3479.8681]
        {'alt1', 400, 'alt2', 600, 'apo2', 39700, 'elements', [51.6 30 0 63.4 100 270], 'tof_max_h', 30}, ...
        [58.815692017167976 0 210.591986530302], [4.442609368 8948.9431]};
for q = 1:size(pyEl, 1)
    m = odst_mission(pyEl{q, 1}{:});
    opts = odst_evaluate(m); opts = opts([opts.applicable]);
    got = [[opts.dv].', [opts.tof].'];
    report(max(abs([m.di_deg m.w1_deg m.dw_deg] - pyEl{q, 2})) < 1e-9 && isequal(size(got), size(pyEl{q, 3})) && ...
           max(abs(got(:, 1) - pyEl{q, 3}(:, 1))) < 1e-6 && max(abs(got(:, 2) - pyEl{q, 3}(:, 2))) < 1, ...
           sprintf('[%s]: di %.6f, w1 %.6f, dw %.6f deg and %d options as Python (dv to 1e-6 km/s, TOF to 1 s)', ...
                   strtrim(sprintf('%g ', m.elements)), m.di_deg, m.w1_deg, m.dw_deg, numel(opts)));
end
% Every plan, turned into EarthMJ2000Eq and propagated, arrives on the target given by its own elements
elm = {{'alt1', 400, 'alt2', 35786, 'elements', [51.6 120 0 0 0 0]}
       {'alt1', 400, 'alt2', 35786, 'elements', [28.5 0 0 0 0 0], 'tof_max_h', 3}
       {'alt1', 500, 'alt2', 500, 'elements', [51.6 0 0 51.6 10 0]}
       {'alt1', 250, 'apo1', 35786, 'alt2', 35786, 'elements', [6 0 178 0 0 0]}
       {'alt1', 400, 'alt2', 600, 'apo2', 39700, 'elements', [51.6 30 0 63.4 100 270], 'tof_max_h', 30}
       {'alt1', 300, 'apo1', 12000, 'alt2', 800, 'apo2', 20000, 'elements', [28.5 10 30 40 70 120], 'tof_max_h', 20}};
for q = 1:numel(elm)
    m = odst_mission(elm{q}{:});
    e = m.elements*pi/180;
    [P2, ~, W2] = odst_perifocal(e(4), e(5), e(6));
    rp = RE + m.alt2; ra = RE + max(m.alt2, m.apo2); p = 2*rp*ra/(rp + ra);
    opts = odst_evaluate(m);
    for i = 1:numel(opts)
        o = opts(i);
        if isempty(o.vplan) || ~o.applicable, continue; end
        [y0i, plani] = odst_inertial_plan(m, o);
        [~, ~, yf, dvSim] = odst_simulate_vec(y0i, plani, 3600);
        r = yf(1:3); v = yf(4:6); h = cross(r, v);
        dh = norm(h - sqrt(C.MU*p)*W2)/sqrt(C.MU*p);
        de = norm(cross(v, h)/C.MU - r/norm(r) - (ra - rp)/(ra + rp)*P2);
        report(dh < 1e-6 && de < 1e-6 && abs(dvSim - o.dv) < 1e-9, sprintf(['[%-22s] %-38s dv=%.6f  arrival vs the ' ...
               'target''s elements |dh|/h=%.1e |de|=%.1e'], strtrim(sprintf('%g ', m.elements)), o.name, o.dv, dh, de));
    end
end

if nFail == 0
    fprintf('\nALL CHECKS PASSED\n');
else
    fprintf('\n%d CHECK(S) FAILED\n', nFail);
end

    function check(name, got, ref, tol)
        report(abs(got - ref) <= tol, sprintf('%s: got %.5f, ref %.5f (tol %g)', name, got, ref, tol));
    end
    function report(ok, msg)
        if ok, tag = 'PASS'; else, tag = 'FAIL'; nFail = nFail + 1; end
        fprintf('[%s] %s\n', tag, msg);
    end
end

function d = evalDv(opts)
% dv of every option, NaN where not applicable
d = [opts.dv]; d(~[opts.applicable]) = NaN;
d(isnan(d)) = 0;
end

function [dh, de] = targetMismatch(y, m)
% Relative angular-momentum error and absolute eccentricity-vector error of
% state y against the mission's target orbit (as target_mismatch in odst.py).
C = odst_const();
rp = C.RE + m.alt2; ra = C.RE + max(m.alt2, m.apo2);
[~, ~, P, Q] = odst_orbit_state(rp, ra, m.di_deg*pi/180, m.dw_deg*pi/180, 0);
r = y(1:3); v = y(4:6);
h = cross(r, v); e = cross(v, h)/C.MU - r/norm(r);
ht = sqrt(C.MU*2*rp*ra/(rp + ra))*cross(P, Q);
dh = norm(h - ht)/norm(ht);
de = norm(e - (ra - rp)/(ra + rp)*P);
end
