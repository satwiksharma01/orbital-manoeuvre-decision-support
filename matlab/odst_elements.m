function [a, e, incDeg] = odst_elements(y)
%ODST_ELEMENTS Semi-major axis [km], eccentricity and inclination [deg] of state y.
C = odst_const(); mu = C.MU;
r = y(1:3); v = y(4:6);
r = r(:); v = v(:);
h = cross(r, v);
ev = cross(v, h)/mu - r/norm(r);
a = 1/(2/norm(r) - dot(v, v)/mu);
e = norm(ev);
incDeg = acos(h(3)/norm(h))*180/pi;
end
