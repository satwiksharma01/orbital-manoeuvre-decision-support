function [di, w1, dw, R] = odst_relative_geometry(m)
%ODST_RELATIVE_GEOMETRY Relative geometry of a mission given by orbital elements.
%   [di, w1, dw, R] = odst_relative_geometry(m), with m.elements = [i1 RAAN1 argp1
%   i2 RAAN2 argp2] in deg: the angle between the planes, each periapsis measured
%   from the node line, where the target crosses the start plane going north [deg],
%   and R, whose columns are the relative frame's axes in EarthMJ2000Eq. Coplanar
%   orbits have no node line, so the start periapsis takes its place; a circular
%   orbit's periapsis is put on the node line. Mirrors relative_geometry() in
%   python/odst.py.
e = m.elements*pi/180;
[P1, ~, W1] = odst_perifocal(e(1), e(2), e(3));
[P2, ~, W2] = odst_perifocal(e(4), e(5), e(6));
n = cross(W1, W2);
if norm(n) > 1e-9, x = n/norm(n); else, x = P1; end
y = cross(W1, x); y2 = cross(W2, x);              % 90 deg ahead of +x in each plane
di = atan2(norm(n), dot(W1, W2))*180/pi;
w1 = fromNode(P1, y, isnan(m.apo1) || m.apo1 - m.alt1 < 1e-6);
dw = fromNode(P2, y2, isnan(m.apo2) || m.apo2 - m.alt2 < 1e-6);
R = [x, y, W1];

    function a = fromNode(P, yy, circ)
        a = 0;                                    % round: -1e-16 deg is 0, not 360
        if ~circ, a = mod(round(atan2(dot(P, yy), dot(P, x))*180/pi, 12), 360); end
    end
end
