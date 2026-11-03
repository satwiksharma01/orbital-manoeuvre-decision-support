function [y0, burns] = odst_local_burns(m, o)
%ODST_LOCAL_BURNS With J2: option o's start state and its burns in each burn's local VNB frame.
%   [y0, burns] = odst_local_burns(m, o): y0 in EarthMJ2000Eq; burns, a struct array with
%   t (planned time, s), dv (1 x 3 km/s along velocity, orbit normal and binormal), at ('time':
%   after the planned time; 'apoapsis' or 'periapsis': the burn ends a coast there; 'node': after
%   the planned wait, at the next equator crossing) and signed (the normal component is for an
%   ascending node and flips at a descending one). A burn's direction then follows the plane J2
%   turns, which a fixed EarthMJ2000Eq vector doesn't. Mirrors local_burns() in python/odst.py.
C = odst_const();
burns = struct('t', {}, 'dv', {}, 'at', {}, 'signed', {});
if ~strcmp(o.key, 'nodal_drift')
    [y0, vplan] = odst_inertial_plan(m, o);
    y = y0; t0 = 0;
    for k = 1:size(vplan, 1)
        t = vplan(k, 1); dv = vplan(k, 2:4).';
        if t > t0, y = odst_kepler(y, t - t0); end
        r = y(1:3); v = y(4:6); at = 'time';
        if t > 0 && abs(dot(r, v)) < 1e-6*norm(r)*norm(v)
            if norm(r) > 1/(2/norm(r) - dot(v, v)/C.MU), at = 'apoapsis'; else, at = 'periapsis'; end
        end
        n = cross(r, v)/norm(cross(r, v)); V = v/norm(v);
        burns(end+1) = struct('t', t, 'dv', [dot(V, dv), dot(n, dv), dot(cross(V, n), dv)], 'at', at, 'signed', false); %#ok<AGROW>
        y(4:6) = v + dv; t0 = t;
    end
    return
end
[P, Q] = odst_perifocal(m.elements(1)*pi/180, m.elements(2)*pi/180, 0);
r1 = C.RE + m.alt1;
y0 = [r1*P; sqrt(C.MU/r1)*Q];                  % on the start orbit's ascending node
s = 1; if m.elements(4) < m.elements(1), s = -1; end    % the sign of the inclination change
rp = r1; vp = sqrt(C.MU/r1);
for k = 1:size(o.plan, 1)
    t = o.plan(k, 1); va = o.plan(k, 2); dinc = o.plan(k, 3); r = o.plan(k, 4);
    vb = sqrt(C.MU*(2/r - (2/rp - vp*vp/C.MU)));      % the speed before the burn, on the orbit flown since the last
    if t == 0
        at = 'time';
    elseif abs(r - rp) < 1e-6
        if o.params.di_deg > 0, at = 'node'; else, at = 'time'; end
    elseif r > rp
        at = 'apoapsis';
    else
        at = 'periapsis';
    end
    burns(end+1) = struct('t', t, 'dv', [va*cos(dinc) - vb, s*va*sin(dinc), 0], 'at', at, 'signed', dinc > 0); %#ok<AGROW>
    rp = r; vp = va;
end
end
