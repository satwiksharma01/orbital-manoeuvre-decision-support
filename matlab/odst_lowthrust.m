function [dv, tof] = odst_lowthrust(r1, r2, di, accel_ms2, isp)
%ODST_LOWTHRUST Edelbaum continuous low-thrust transfer between circular orbits.
%   dv = sqrt(v1^2 + v2^2 - 2 v1 v2 cos(pi/2 * di)); di in rad, accel_ms2 the
%   acceleration at departure in m/s^2. Without isp it stays constant:
%   TOF = dv / a. With isp the thrust is constant and the mass falls:
%   TOF = ve / a (1 - exp(-dv / ve)); Edelbaum's dv and yaw law are unchanged, as
%   both depend only on the velocity increment u = integral of F/m dt. Valid for
%   di <= 2 rad (114.6 deg); beyond it dv falls as di grows, so odst_evaluate
%   marks the option not applicable. Mirrors low_thrust() in python/odst.py.
C = odst_const(); mu = C.MU;
v1 = sqrt(mu/r1); v2 = sqrt(mu/r2);
dv = sqrt(max(v1^2 + v2^2 - 2*v1*v2*cos(pi/2*di), 0));
if nargin < 5 || isempty(isp)
    tof = dv/(accel_ms2*1e-3);
else
    ve = isp*C.G0;
    tof = ve/(accel_ms2*1e-3)*(1 - exp(-dv/ve));
end
end
