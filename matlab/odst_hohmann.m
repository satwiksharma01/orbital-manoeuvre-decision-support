function [dv, tof, nburns, plan] = odst_hohmann(r1, r2, di, split)
%ODST_HOHMANN Hohmann transfer between circular orbits with a plane change.
%   [dv, tof, nburns, plan] = odst_hohmann(r1, r2, di, split)
%   r1, r2 [km], di [rad]. split = fraction of di performed at the r1 burn.
%   split = [] (default) applies the classical rule: whole plane change at the
%   burn on the larger radius, where the orbital speed is lowest.
%   plan = [t_burn, speed_after_burn, plane_rotation] per burn (one row each).
C = odst_const(); mu = C.MU;
if nargin < 4 || isempty(split)
    split = double(r2 < r1);
end
if abs(r2 - r1) < 1e-6                 % degenerate: pure plane change at a node
    v = sqrt(mu/r1);
    dv = 2*v*sin(di/2); tof = 0; nburns = 1; plan = [0, v, di];
    return
end
a   = (r1 + r2)/2;
v1  = sqrt(mu/r1);              v2  = sqrt(mu/r2);
vt1 = sqrt(mu*(2/r1 - 1/a));    vt2 = sqrt(mu*(2/r2 - 1/a));
dv1 = odst_combined(v1, vt1, split*di);
dv2 = odst_combined(vt2, v2, (1 - split)*di);
dv  = dv1 + dv2;
tof = pi*sqrt(a^3/mu);
nburns = 2;
plan = [0,   vt1, split*di;
        tof, v2,  (1 - split)*di];
end
