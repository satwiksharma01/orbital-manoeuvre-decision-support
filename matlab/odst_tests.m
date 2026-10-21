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
       odst_mission('alt1', 500, 'apo1', 5000, 'alt2', 800, 'apo2', 12000, 'dw_deg', 60, 'di_deg', 10)};
pyRef = {[1.833344084 0; 1.833344097 20504.3295], [2.476172306 21378.6762; 2.476173769 21646.7968], ...
         [NaN NaN; 1.183670980 9359.7970], [NaN NaN; 1.859013203 3601.2405]};
pyRec = {'apse', 'apse', 'lambert', 'lambert'};
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
        report(abs(o.dv - pyRef{c}(i, 1)) < 1e-6 && abs(o.tof - pyRef{c}(i, 2)) < 1, ...
               sprintf('%s: %-8s dv %.9f vs Python %.9f (%.1e km/s), TOF %.4f vs %.4f s', mt, o.key, ...
                       o.dv, pyRef{c}(i, 1), abs(o.dv - pyRef{c}(i, 1)), o.tof, pyRef{c}(i, 2)));
        [~, ~, yf, dvSim] = odst_simulate_vec(o.y0, o.vplan, 3600);
        [dh, de] = targetMismatch(yf, m);
        report(dh < 1e-6 && de < 1e-6 && abs(dvSim - o.dv) < 1e-9, ...
               sprintf('%s: %-8s propagated arrival |dh|/h = %.1e, |de| = %.1e', mt, o.key, dh, de));
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
