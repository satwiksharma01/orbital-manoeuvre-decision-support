function [T, Y, yf, dvUsed] = odst_simulate_vec(y0, plan, coastAfter)
%ODST_SIMULATE_VEC Verify a vector burn plan by two-body numerical propagation.
%   [T, Y, yf, dvUsed] = odst_simulate_vec(y0, plan, coastAfter) starts from
%   state y0 = [r; v] and adds each burn of plan = [t, dvx, dvy, dvz] at its
%   time. Used for every impulsive option's vector plan (odst_evaluate);
%   mirrors simulate_vec() in python/odst.py.
if nargin < 3, coastAfter = 0; end
C = odst_const(); mu = C.MU;
f = @(~, y) [y(4:6); -mu*y(1:3)/norm(y(1:3))^3];
ode = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);
y = y0(:); t0 = 0; T = []; Y = []; dvUsed = 0;
nb = size(plan, 1);
if nb > 0, tEnd = plan(end, 1) + coastAfter; else, tEnd = coastAfter; end
times = [plan(:, 1); tEnd];
for i = 1:nb + 1
    tb = times(i);
    if tb > t0
        [tt, yy] = ode45(f, [t0, tb], y, ode);
        T = [T; tt]; Y = [Y; yy]; %#ok<AGROW>
        y = yy(end, :).';
    end
    if i <= nb
        y(4:6) = y(4:6) + plan(i, 2:4).';
        dvUsed = dvUsed + norm(plan(i, 2:4));
    end
    t0 = tb;
end
if isempty(T), T = 0; Y = y.'; end
yf = y;
end
