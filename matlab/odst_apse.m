function [res, why] = odst_apse(rp1, ra1, rp2, ra2, di, dw)
%ODST_APSE Apse-to-apse transfer with the optimal plane-change split.
%   [res, why] = odst_apse(rp1, ra1, rp2, ra2, di, dw): two impulses, from an
%   apsis of the start orbit to the opposite direction on the target orbit,
%   with the plane change di split between the burns by fminbnd. Burns below
%   1e-9 km/s are dropped, so GTO -> GEO is one apogee burn. res = [] with the
%   reason in why when the apse lines don't allow it. Radii in km, angles in
%   rad (geometry as in odst_orbit_state). Mirrors apse_transfer() in
%   python/odst.py.
C = odst_const(); mu = C.MU;
res = []; why = '';
circ1 = ra1 - rp1 < 1e-6; circ2 = ra2 - rp2 < 1e-6;
onNode = @(a) abs(sin(a)) < 1e-9;
if ~(circ1 || circ2 || onNode(dw))
    why = 'apse lines not aligned (dw must be 0 or 180 deg)';
    return
end
if circ1 && ~circ2, dirs = [dw, dw + pi]; else, dirs = [0, pi]; end
if di > 0 && ~(onNode(dirs(1)) && onNode(dirs(2)))
    why = 'a plane change needs both burns on the node line (dw must be 0 or 180 deg)';
    return
end
vv = @(r, a) sqrt(mu*(2/r - 1/a));
a1 = (rp1 + ra1)/2; a2 = (rp2 + ra2)/2;
best = Inf;
for phi = dirs
    th2 = phi + pi - dw;
    [~, ~, ~, ~, rA] = odst_orbit_state(rp1, ra1, 0, 0, phi);
    [~, ~, ~, ~, rB] = odst_orbit_state(rp2, ra2, di, dw, th2);
    at = (rA + rB)/2;
    v1 = vv(rA, a1); vtA = vv(rA, at); vtB = vv(rB, at); v2 = vv(rB, a2);
    cost = @(s) odst_combined(v1, vtA, s*di) + odst_combined(vtB, v2, (1 - s)*di);
    cands = [0 1];
    if di > 0
        cands(3) = fminbnd(cost, 0, 1, optimset('TolX', 1e-10));
    end
    [cmin, j] = min(arrayfun(cost, cands));
    if cmin < best - 1e-12
        best = cmin;
        B = struct('phi', phi, 'th2', th2, 'rA', rA, 'rB', rB, 'at', at, 'v1', v1, ...
                   'vtA', vtA, 'vtB', vtB, 'v2', v2, 's', cands(j));
    end
end
s = B.s; phi = B.phi;
dv1 = odst_combined(B.v1, B.vtA, s*di);
dv2 = odst_combined(B.vtB, B.v2, (1 - s)*di);
Rt = [1 0 0; 0 cos(s*di) -sin(s*di); 0 sin(s*di) cos(s*di)];   % transfer plane
along = @(ang) Rt*[-sin(ang); cos(ang); 0];
[r1v, v1v] = odst_orbit_state(rp1, ra1, 0, 0, phi);
rBv = Rt*[B.rB*cos(phi + pi); B.rB*sin(phi + pi); 0];
vB = B.vtB*along(phi + pi);
dv1v = B.vtA*along(phi) - v1v;
[~, v2t] = odst_orbit_state(rp2, ra2, di, dw, B.th2);
dv2v = v2t - vB;
tof = pi*sqrt(B.at^3/mu);
keep1 = dv1 > 1e-9; keep2 = dv2 > 1e-9;
if keep1 && keep2
    y0 = [r1v; v1v]; plan = [0, dv1v.'; tof, dv2v.']; legs = [B.rA, dv1; B.rB, dv2];
elseif keep2                       % already on the transfer orbit: one burn on arrival
    y0 = [rBv; vB]; plan = [0, dv2v.']; legs = [B.rB, dv2]; tof = 0;
elseif keep1                       % the first burn lands on the target orbit
    y0 = [r1v; v1v]; plan = [0, dv1v.']; legs = [B.rA, dv1]; tof = 0;
else
    y0 = [r1v; v1v]; plan = zeros(0, 4); legs = zeros(0, 2); tof = 0;
end
res = struct('dv', dv1*keep1 + dv2*keep2, 'tof', tof, 'burns', size(plan, 1), 'y0', y0, ...
             'plan', plan, 'legs', legs, 'split_at_dep', s, 'depart_deg', mod(phi*180/pi, 360), ...
             'rA', B.rA, 'rB', B.rB);
end
