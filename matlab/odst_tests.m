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
