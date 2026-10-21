function res = odst_lambert_opt(rp1, ra1, rp2, ra2, di, dw, tofmax, nGrid)
%ODST_LAMBERT_OPT Minimum-dv two-impulse Lambert transfer between two orbits.
%   res = odst_lambert_opt(rp1, ra1, rp2, ra2, di, dw, tofmax, nGrid) searches
%   departure true anomaly, arrival true anomaly and transfer time (60 s to
%   min(tofmax, the longer orbital period)): an nGrid x nGrid grid (default
%   16) with a golden-section search on time per cell, then Nelder-Mead from
%   the best cell. Transfers within 0.5 deg of 0 or 180 deg are skipped, as
%   the transfer plane is undefined there. res = [] if nothing is found.
%   Same algorithm, step for step, as lambert_optimal() in python/odst.py,
%   so the engines agree to rounding.
if nargin < 8, nGrid = 16; end
C = odst_const(); mu = C.MU;
tmin = 60;
per = @(rp, ra) 2*pi*sqrt(((rp + ra)/2)^3/mu);
tmax = min(tofmax, max(per(rp1, ra1), per(rp2, ra2)));
res = [];
if tmax <= tmin, return; end
gap = sin(0.5*pi/180);
cost = @(x) solve(x);
g = (sqrt(5) - 1)/2;
fBest = Inf; xBest = [];
for i = 0:nGrid - 1
    for j = 0:nGrid - 1
        th1 = 2*pi*i/nGrid; th2 = 2*pi*j/nGrid;
        a = 0; b = 1;
        c = b - g*(b - a); d = a + g*(b - a);
        fc = cost([th1, th2, c]); fd = cost([th1, th2, d]);
        for k = 1:25
            if fc < fd
                b = d; d = c; fd = fc;
                c = b - g*(b - a);
                fc = cost([th1, th2, c]);
            else
                a = c; c = d; fc = fd;
                d = a + g*(b - a);
                fd = cost([th1, th2, d]);
            end
        end
        u = (a + b)/2;
        fu = cost([th1, th2, u]);
        if fu < fBest
            fBest = fu; xBest = [th1, th2, u];
        end
    end
end
if ~isfinite(fBest), return; end
step = 2*pi/nGrid;
[x, fx] = nelderMead(cost, xBest, [step/2, step/2, 0.05]);
[~, r1v, v1v, r2v, dv1v, dv2v] = solve(x);
T = tmin + x(3)*(tmax - tmin);
res = struct('dv', fx, 'tof', T, 'burns', 2, 'y0', [r1v; v1v], 'plan', [0, dv1v.'; T, dv2v.'], ...
             'legs', [norm(r1v), norm(dv1v); norm(r2v), norm(dv2v)], ...
             'th1_deg', mod(x(1)*180/pi, 360), 'th2_deg', mod(x(2)*180/pi, 360), ...
             'time_limited', x(3) > 1 - 1e-6);

    function [f, r1v, v1v, r2v, dv1v, dv2v] = solve(x)
        f = Inf; r1v = []; v1v = []; r2v = []; dv1v = []; dv2v = [];
        if x(3) < 0 || x(3) > 1, return; end
        [r1v, v1v] = odst_orbit_state(rp1, ra1, 0, 0, x(1));
        [r2v, v2v] = odst_orbit_state(rp2, ra2, di, dw, x(2));
        if norm(cross(r1v, r2v)) < gap*norm(r1v)*norm(r2v), return; end
        [w1, w2] = odst_lambert(r1v, r2v, tmin + x(3)*(tmax - tmin));
        if isempty(w1), return; end
        dv1v = w1 - v1v; dv2v = v2v - w2;
        f = norm(dv1v) + norm(dv2v);
    end
end

function [xb, fb] = nelderMead(f, x0, steps)
% Plain Nelder-Mead (reflect 1, expand 2, contract 0.5, shrink 0.5), the same
% steps as _nelder_mead() in python/odst.py.
n = numel(x0); iters = 400; tol = 1e-12;
pts = repmat(x0, n + 1, 1);
for i = 1:n, pts(i + 1, i) = pts(i + 1, i) + steps(i); end
vals = zeros(n + 1, 1);
for i = 1:n + 1, vals(i) = f(pts(i, :)); end
for it = 1:iters
    [vals, o] = sort(vals); pts = pts(o, :);
    if abs(vals(end) - vals(1)) < tol, break; end
    c = sum(pts(1:n, :), 1)/n;
    xr = c + (c - pts(end, :)); fr = f(xr);
    if fr < vals(1)
        xe = c + 2*(c - pts(end, :)); fe = f(xe);
        if fe < fr
            pts(end, :) = xe; vals(end) = fe;
        else
            pts(end, :) = xr; vals(end) = fr;
        end
    elseif fr < vals(end - 1)
        pts(end, :) = xr; vals(end) = fr;
    else
        xc = c + 0.5*(pts(end, :) - c); fc = f(xc);
        if fc < vals(end)
            pts(end, :) = xc; vals(end) = fc;
        else
            for i = 2:n + 1
                pts(i, :) = pts(1, :) + 0.5*(pts(i, :) - pts(1, :));
                vals(i) = f(pts(i, :));
            end
        end
    end
end
[fb, i] = min(vals);
xb = pts(i, :);
end
