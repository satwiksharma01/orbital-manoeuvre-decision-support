function [v1, v2] = odst_lambert(r1v, r2v, dt)
%ODST_LAMBERT Single-revolution prograde Lambert problem.
%   [v1, v2] = odst_lambert(r1v, r2v, dt) returns the velocities [km/s] at r1v
%   and r2v [km] of the transfer that takes dt [s], with its angular momentum
%   along +z. Universal variables (Curtis Algorithm 5.2) with bisection on z.
%   v1 = v2 = [] when there is no solution. Mirrors lambert() in python/odst.py.
C = odst_const(); mu = C.MU;
r1v = r1v(:); r2v = r2v(:);
r1 = norm(r1v); r2 = norm(r2v);
v1 = []; v2 = [];
th = acos(max(-1, min(1, dot(r1v, r2v)/(r1*r2))));
if r1v(1)*r2v(2) - r1v(2)*r2v(1) < 0
    th = 2*pi - th;
end
if 1 - cos(th) < 1e-14, return; end
A = sin(th)*sqrt(r1*r2/(1 - cos(th)));
if abs(A) < 1e-12, return; end

lo = -4*pi; hi = 4*pi^2 - 1e-4;            % C(z) -> 0 at 4 pi^2
while tofOf(lo, r1, r2, A, mu) > dt        % very fast: needs a more hyperbolic bracket
    lo = 2*lo;
    if lo < -1e5, return; end
end
if tofOf(hi, r1, r2, A, mu) < dt, return; end
for it = 1:200
    z = 0.5*(lo + hi);
    if tofOf(z, r1, r2, A, mu) < dt
        lo = z;
    else
        hi = z;
    end
    if hi - lo < 1e-13*max(1, abs(z)), break; end
end
y = yOf(0.5*(lo + hi), r1, r2, A);
f = 1 - y/r1; g = A*sqrt(y/mu); gdot = 1 - y/r2;
v1 = (r2v - f*r1v)/g;
v2 = (gdot*r2v - r1v)/g;
end

function [y, Cz, Sz] = yOf(z, r1, r2, A)
[Cz, Sz] = stumpff(z);
y = r1 + r2 + A*(z*Sz - 1)/sqrt(Cz);
end

function t = tofOf(z, r1, r2, A, mu)
[y, Cz, Sz] = yOf(z, r1, r2, A);
if y < 0
    t = -Inf;
    return
end
x = sqrt(y/Cz);
t = (x^3*Sz + A*sqrt(y))/sqrt(mu);
end

function [Cz, Sz] = stumpff(z)
if z > 1e-8
    s = sqrt(z);
    Cz = (1 - cos(s))/z; Sz = (s - sin(s))/s^3;
elseif z < -1e-8
    s = sqrt(-z);
    Cz = (cosh(s) - 1)/-z; Sz = (sinh(s) - s)/s^3;
else
    Cz = 0.5 - z/24; Sz = 1/6 - z/120;
end
end
