function [dv, tof, nburns, plan, info] = odst_onetangent(r1, r2, k, di)
%ODST_ONETANGENT One-tangent-burn (fast) transfer between circular orbits.
%   Raising (r2 > r1): tangential departure at periapsis r1, apoapsis k*r2.
%   Lowering (r2 < r1): tangential departure at apoapsis r1, periapsis r2/k.
%   k >= 1 (k = 1 is the Hohmann transfer). The arrival burn is not tangential
%   (flight-path angle phi). The plane change is combined with the burn on the
%   larger radius. info: e, true anomaly at arrival nu [deg], phi [deg].
C = odst_const(); mu = C.MU;
raising = r2 > r1;
if raising
    rp = r1;  ra = k*r2;
else
    ra = r1;  rp = r2/k;
end
a = (rp + ra)/2;
e = (ra - rp)/(ra + rp);
info = struct('e', e, 'nu_deg', NaN, 'fpa_deg', NaN, 'k', k);
if e < 1e-12
    [dv, tof, nburns, plan] = odst_hohmann(r1, r2, di);
    return
end
p  = a*(1 - e^2);
v1 = sqrt(mu/r1);              v2 = sqrt(mu/r2);
vt1 = sqrt(mu*(2/r1 - 1/a));   vt2 = sqrt(mu*(2/r2 - 1/a));
nu  = acos(min(1, max(-1, (p/r2 - 1)/e)));          % 0..pi
phi = atan2(e*sin(nu), 1 + e*cos(nu));               % flight-path angle
E   = 2*atan(sqrt((1 - e)/(1 + e))*tan(nu/2));
tPeri = sqrt(a^3/mu)*(E - e*sin(E));                 % time from periapsis
if raising
    tof = tPeri;
    dv1 = abs(vt1 - v1);
    dv2 = sqrt(max(vt2^2 + v2^2 - 2*vt2*v2*cos(phi)*cos(di), 0));
    plan = [0, vt1, 0, r1, dv1; tof, v2, di, r2, dv2];
else
    tof = pi*sqrt(a^3/mu) - tPeri;
    dv1 = odst_combined(v1, vt1, di);
    dv2 = sqrt(max(vt2^2 + v2^2 - 2*vt2*v2*cos(phi), 0));
    plan = [0, vt1, di, r1, dv1; tof, v2, 0, r2, dv2];
end
dv = dv1 + dv2;
nburns = 2;
info.nu_deg = nu*180/pi;
info.fpa_deg = phi*180/pi;
end
