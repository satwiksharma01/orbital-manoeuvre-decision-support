function [P, Q, W] = odst_perifocal(inc, raan, argp)
%ODST_PERIFOCAL Unit vectors to periapsis (P), 90 deg ahead of it (Q) and along the
%   orbit normal (W) in EarthMJ2000Eq, from inclination, RAAN and argument of
%   periapsis [rad]. Mirrors perifocal() in python/odst.py.
ci = cos(inc); si = sin(inc); cO = cos(raan); sO = sin(raan); cw = cos(argp); sw = sin(argp);
P = [cO*cw - sO*sw*ci; sO*cw + cO*sw*ci; sw*si];
Q = [-cO*sw - sO*cw*ci; -sO*sw + cO*cw*ci; cw*si];
W = [sO*si; -cO*si; ci];
end
