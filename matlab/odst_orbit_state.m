function [r, v, P, Q, rr] = odst_orbit_state(rp, ra, di, dw, th)
%ODST_ORBIT_STATE Position [km] and velocity [km/s] at true anomaly th [rad]
%   on the orbit with periapsis and apoapsis radii rp, ra [km]. Its plane is
%   turned by di about +x, and its periapsis lies dw from +x in that plane
%   (the start orbit is di = dw = 0). P, Q: unit vectors to periapsis and 90
%   deg ahead of it; rr: the radius, exact at the apsides.
%   Mirrors orbit_basis/orbit_radius/orbit_state in python/odst.py.
C = odst_const();
x = [1; 0; 0]; y = [0; cos(di); sin(di)];
P = cos(dw)*x + sin(dw)*y;
Q = -sin(dw)*x + cos(dw)*y;
e = (ra - rp)/(ra + rp); p = 2*rp*ra/(rp + ra);
if abs(cos(th) - 1) < 1e-12
    rr = rp;
elseif abs(cos(th) + 1) < 1e-12
    rr = ra;
else
    rr = p/(1 + e*cos(th));
end
r = rr*(cos(th)*P + sin(th)*Q);
v = sqrt(C.MU/p)*(-sin(th)*P + (e + cos(th))*Q);
end
