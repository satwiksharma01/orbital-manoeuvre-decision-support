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
