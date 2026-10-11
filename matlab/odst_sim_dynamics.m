function xdot = odst_sim_dynamics(u)
%ODST_SIM_DYNAMICS Two-body dynamics + finite-burn thrust for the Simulink model.
%   u = [t; r(3); v(3)].  Burn plan read from base workspace variable
%   odst_burnplan (fields t, vafter, dinc, tau, mu). Each impulsive burn k is
%   replaced by a constant acceleration dv_k/tau over [t_k - tau/2, t_k + tau/2].
%   The burn vector is fixed in inertial space when the burn starts. It is
%   targeted at the burn CENTRE: the current state is coasted to t_k (RK4), the
%   post-burn velocity there is built (speed vafter, perpendicular to r, in the
%   plane rotated by dinc about r) and dv_k = v_target - v_centre.
%   Call "clear odst_sim_dynamics" before each run.
persistent P acc started
if isempty(P)
    P = evalin('base', 'odst_burnplan');
    acc = zeros(numel(P.t), 3);
    started = false(numel(P.t), 1);
end
t = u(1);
r = u(2:4); r = r(:);
v = u(5:7); v = v(:);
a = -P.mu*r/norm(r)^3;
for k = 1:numel(P.t)
    if t >= P.t(k) - P.tau/2 && t < P.t(k) + P.tau/2
        if ~started(k)
            [rc, vc] = coast(r, v, P.t(k) - t, P.mu);
            h = cross(rc, vc); nh = h/norm(h);
            nA = rodrigues(nh, rc, P.dinc(k));
            nB = rodrigues(nh, rc, -P.dinc(k));
            if nA(3) <= nB(3), nNew = nA; else, nNew = nB; end
            vTarget = P.vafter(k)*cross(nNew, rc/norm(rc));
            acc(k, :) = ((vTarget - vc)/P.tau).';
            started(k) = true;
        end
        a = a + acc(k, :).';
    end
end
xdot = [v; a];
end

function v = rodrigues(v, axis, ang)
k = axis/norm(axis);
v = v*cos(ang) + cross(k, v)*sin(ang) + k*dot(k, v)*(1 - cos(ang));
end

function [r, v] = coast(r, v, dt, mu)
% Fixed-step RK4 two-body coast over dt (50 steps).
n = 50; h = dt/n; y = [r; v];
g = @(y) [y(4:6); -mu*y(1:3)/norm(y(1:3))^3];
for i = 1:n
    k1 = g(y); k2 = g(y + h/2*k1); k3 = g(y + h/2*k2); k4 = g(y + h*k3);
    y = y + h/6*(k1 + 2*k2 + 2*k3 + k4);
end
r = y(1:3); v = y(4:6);
end
