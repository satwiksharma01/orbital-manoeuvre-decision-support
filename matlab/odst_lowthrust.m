function [dv, tof] = odst_lowthrust(r1, r2, di, accel_ms2)
%ODST_LOWTHRUST Edelbaum continuous low-thrust transfer between circular orbits.
%   dv = sqrt(v1^2 + v2^2 - 2 v1 v2 cos(pi/2 * di)),  TOF = dv / (F/m)
%   (constant-acceleration approximation; di in rad, accel in m/s^2). Valid for
%   di <= 2 rad (114.6 deg); beyond it dv falls as di grows, so odst_evaluate
%   marks the option not applicable.
C = odst_const(); mu = C.MU;
v1 = sqrt(mu/r1); v2 = sqrt(mu/r2);
dv = sqrt(max(v1^2 + v2^2 - 2*v1*v2*cos(pi/2*di), 0));
tof = dv/(accel_ms2*1e-3);
end
