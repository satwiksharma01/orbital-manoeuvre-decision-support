function y = odst_kepler(y, dt)
%ODST_KEPLER Two-body state dt seconds after state y (6x1, km and km/s).
%   The universal-variable Kepler equation is solved by Newton, then the state follows
%   from the Lagrange f and g functions (Curtis, Algorithms 3.3 and 3.4), as kepler in odst.py.
C = odst_const(); mu = C.MU; sq = sqrt(mu);
r0v = y(1:3); v0v = y(4:6);
r0 = norm(r0v); vr0 = dot(r0v, v0v)/r0; alpha = 2/r0 - dot(v0v, v0v)/mu;
x = sq*abs(alpha)*dt;
for it = 1:100
    z = alpha*x^2; [Cz, Sz] = stumpff(z);
    F = r0*vr0/sq*x^2*Cz + (1 - alpha*r0)*x^3*Sz + r0*x - sq*dt;
    dx = F/(r0*vr0/sq*x*(1 - z*Sz) + (1 - alpha*r0)*x^2*Cz + r0);
    x = x - dx;
    if abs(dx) < 1e-12*max(1, abs(x)), break; end
    if it == 100, error('odst_kepler:noConvergence', 'Newton did not converge.'); end
end
z = alpha*x^2; [Cz, Sz] = stumpff(z);
rv = (1 - x^2/r0*Cz)*r0v + (dt - x^3*Sz/sq)*v0v;
r = norm(rv);
y = [rv; sq/(r*r0)*x*(z*Sz - 1)*r0v + (1 - x^2/r*Cz)*v0v];
end

function [Cz, Sz] = stumpff(z)
if z > 1e-8
    s = sqrt(z); Cz = (1 - cos(s))/z; Sz = (s - sin(s))/s^3;
elseif z < -1e-8
    s = sqrt(-z); Cz = (cosh(s) - 1)/(-z); Sz = (sinh(s) - s)/s^3;
else
    Cz = 0.5 - z/24; Sz = 1/6 - z/120;
end
end
