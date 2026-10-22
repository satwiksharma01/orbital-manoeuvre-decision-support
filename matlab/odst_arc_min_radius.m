function rmin = odst_arc_min_radius(r1v, v1v, r2v)
%ODST_ARC_MIN_RADIUS Lowest radius [km] on the conic from r1v (velocity v1v)
%   forward to r2v: the transfer's periapsis radius if periapsis lies between
%   them, otherwise the lower end point. Mirrors arc_min_radius() in
%   python/odst.py; odst_lambert_opt rejects arcs below 100 km altitude.
C = odst_const();
r1v = r1v(:); v1v = v1v(:); r2v = r2v(:);
h = cross(r1v, v1v);
ev = cross(v1v, h)/C.MU - r1v/norm(r1v);
e = norm(ev);
rmin = min(norm(r1v), norm(r2v));
if e < 1e-10, return; end
q = cross(h, ev);                         % 90 deg ahead of periapsis, in the direction of motion
nu = @(r) mod(atan2(dot(q, r)/norm(q), dot(ev, r)/e), 2*pi);
if nu(r2v) < nu(r1v)                      % moving forward from r1 to r2 wraps through periapsis
    rmin = min(rmin, dot(h, h)/C.MU/(1 + e));
end
end
