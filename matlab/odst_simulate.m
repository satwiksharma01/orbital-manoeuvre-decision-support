function [T, Y, yf, dvUsed] = odst_simulate(r1, plan, coastAfter)
%ODST_SIMULATE Verify an impulsive burn plan by two-body numerical propagation.
%   [T, Y, yf, dvUsed] = odst_simulate(r1, plan, coastAfter)
%   Starts on a circular equatorial orbit of radius r1 at (r1, 0, 0) and applies
%   each burn of plan = [t_burn, speed_after, plane_rotation]. Every burn in the
%   ODST plans leaves the velocity perpendicular to r (apsis or circularisation),
%   so a burn is fully defined by the speed after it and a rotation of the
%   orbit plane about the radius vector (the node line).
%   T, Y: time history, yf: final state, dvUsed: sum of applied impulses [km/s].
if nargin < 3, coastAfter = 0; end
C = odst_const(); mu = C.MU;
f = @(~, y) [y(4:6); -mu*y(1:3)/norm(y(1:3))^3];
ode = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);
y = [r1; 0; 0; 0; sqrt(mu/r1); 0];
t0 = 0; T = []; Y = []; dvUsed = 0;
for i = 1:size(plan, 1)
    tb = plan(i, 1);
    if tb > t0
        [tt, yy] = ode45(f, [t0, tb], y, ode);
        T = [T; tt]; Y = [Y; yy]; %#ok<AGROW>
        y = yy(end, :).';
    end
    r = y(1:3); v = y(4:6);
    h = cross(r, v); nh = h/norm(h);
    % rotate about the node line in the sense that increases inclination
    % w.r.t. the initial plane (keeps split plane changes additive)
    nA = rodrigues(nh, r, plan(i, 3));
    nB = rodrigues(nh, r, -plan(i, 3));
    if nA(3) <= nB(3), nNew = nA; else, nNew = nB; end
    vNew = plan(i, 2)*cross(nNew, r/norm(r));
    dvUsed = dvUsed + norm(vNew - v);
    y = [r; vNew];
    t0 = tb;
end
if coastAfter > 0
    [tt, yy] = ode45(f, [t0, t0 + coastAfter], y, ode);
    T = [T; tt]; Y = [Y; yy];
    y = yy(end, :).';
end
if isempty(T), T = t0; Y = y.'; end
yf = y;
end

function v = rodrigues(v, axis, ang)
k = axis/norm(axis);
v = v*cos(ang) + cross(k, v)*sin(ang) + k*dot(k, v)*(1 - cos(ang));
end
