function [opts, rec, advice] = odst_evaluate(m)
%ODST_EVALUATE Evaluate all candidate manoeuvres for mission m and rank them.
%   [opts, rec, advice] = odst_evaluate(odst_mission(...))
%
%   opts   struct array, one element per candidate manoeuvre:
%          key, name, dv [km/s], tof [s], burns, isp [s], complexity,
%          applicable, note, plan (burn plan), params, prop_frac (mp/m0),
%          prop_kg (propellant mass, with the spacecraft model),
%          feasible, reasons, chance (share of the Monte Carlo runs that fit the
%          limits, with m.uncertainty), pareto (on the propellant-time Pareto front),
%          j2_dv (with m.j2: dv minus the two-body plan's), m_eff (with m.j2: the
%          mission the option was planned for, its target turned as seen on arrival),
%          score_ws, score_topsis
%   rec    recommended option (highest weighted-sum score among feasible), or []
%   advice cellstr of rule-of-thumb statements (rule-based decision matrix)
%
%   Pipeline:  1) analytical models  2) numerical optimisation of free
%   parameters  3) rule-based feasibility filter  4) multi-criteria scoring.
%   If either orbit is elliptical (apo1/apo2 set), the apse-to-apse and
%   Lambert transfers are sized instead of the five circular-orbit models.
%   With orbital elements and m.j2, every plan allows for J2 turning the planes,
%   and circular missions get the nodal drift option too.
C = odst_const();
r1 = C.RE + m.alt1;
r2 = C.RE + m.alt2;
di = m.di_deg*pi/180;
ra1 = r1; ra2 = r2;
if ~isnan(m.apo1), ra1 = C.RE + m.apo1; end
if ~isnan(m.apo2), ra2 = C.RE + m.apo2; end
ell = ra1 - r1 > 1e-6 || ra2 - r2 > 1e-6;
if ell
    if ra1 < r1 || ra2 < r2
        error('odst_evaluate:apoapsis', 'Apoapsis altitude below periapsis altitude.');
    end
    c = ellipticalOptions(m, C, r1, ra1, r2, ra2, di, '');
else
    c = circularOptions(m, C, r1, r2, di, '');
end

opts = [c{:}];
if m.j2 && ~isempty(m.elements)
    for i = 1:numel(opts), opts(i) = j2Plan(m, C, opts(i)); end
    if ~ell, opts(end+1) = nodalDrift(m, C); end
end
for i = 1:numel(opts)
    if ~isempty(opts(i).plan)
        opts(i).legs = opts(i).plan(:, 4:5);
        if ~strcmp(opts(i).key, 'nodal_drift')   % J2 turns its plane, which no two-body plan holds: it flies in VNB
            mm = m; if ~isempty(opts(i).m_eff), mm = opts(i).m_eff; end
            [opts(i).y0, opts(i).vplan] = circularVecPlan(C.MU, r1, opts(i).plan, mm.di_deg*pi/180);
        end
    end
    if ~isempty(opts(i).legs)
        opts(i).dv_loss = odst_gravity_loss(opts(i).legs, opts(i).isp, m.tw_chem);
    end
    % Tsiolkovsky; the propellant also pays for the finite-burn loss, as the budget check does
    opts(i).prop_frac = 1 - exp(-(opts(i).dv + opts(i).dv_loss)/(opts(i).isp*C.G0));
    opts(i).prop_kg = m.mass0_kg*opts(i).prop_frac;            % NaN without the spacecraft model
end
opts = applyRules(m, opts);
opts = paretoFront(opts);
opts = scoreOptions(m, opts);

rec = [];
f = find([opts.feasible]);
if ~isempty(f)
    [~, j] = max([opts(f).score_ws]);
    rec = opts(f(j));
end
if ell
    advice = ellipticalAdvice(opts(1));
else
    advice = ruleAdvice(m, r1, r2, di);
end
end

% =========================================================================
function c = ellipticalOptions(m, C, r1, ra1, r2, ra2, di, only)
% Hohmann, bi-elliptic, one-tangent and Edelbaum assume circular orbits, so an
% elliptical mission is sized with the apse-to-apse and Lambert models only.
% only: one model's key, or '' for both.
dw = m.dw_deg*pi/180;
c = {};
om1 = m.w1_deg*pi/180;
if isempty(only) || strcmp(only, 'apse')
    [res, why] = odst_apse(r1, ra1, r2, ra2, di, dw, om1);
    name = 'Apse-to-apse (optimal di split)';
    if isempty(res)
        o = newOpt('apse', name, NaN, NaN, 0, m.isp_chem, C.CPLX.apse, [], struct());
        o.applicable = false; o.note = why;
    else
        o = newOpt('apse', name, res.dv, res.tof, res.burns, m.isp_chem, C.CPLX.apse, [], ...
                   struct('split_at_dep', res.split_at_dep, 'depart_deg', res.depart_deg, ...
                          'rA', res.rA, 'rB', res.rB));
        o.y0 = res.y0; o.vplan = res.plan; o.legs = res.legs;
    end
    c{end+1} = o;
end
if ~isempty(only) && ~strcmp(only, 'lambert'), return; end
res = odst_lambert_opt(r1, ra1, r2, ra2, di, dw, m.tof_max_h*3600, [], om1);
name = 'Lambert (minimum dv)';
if isempty(res)
    o = newOpt('lambert', name, NaN, NaN, 2, m.isp_chem, C.CPLX.lambert, [], struct());
    o.applicable = false; o.note = 'no single-revolution Lambert transfer within the time limit';
else
    o = newOpt('lambert', name, res.dv, res.tof, 2, m.isp_chem, C.CPLX.lambert, [], ...
               struct('th1_deg', res.th1_deg, 'th2_deg', res.th2_deg, 'time_limited', res.time_limited));
    o.y0 = res.y0; o.vplan = res.plan; o.legs = res.legs;
    if res.time_limited, o.note = 'time limit sets the transfer time'; end
end
c{end+1} = o;
end

function adv = ellipticalAdvice(apse)
adv = {['Elliptical orbits: Hohmann, bi-elliptic, one-tangent and Edelbaum low thrust assume ' ...
        'circular orbits, so the apse-to-apse and Lambert transfers are sized instead.']};
if apse.applicable
    adv{end+1} = 'Apse lines aligned: a transfer between apsides is available, usually the cheapest.';
else
    adv{end+1} = sprintf('No apse-to-apse transfer: %s. Only a Lambert transfer is offered.', apse.note);
end
end

function c = circularOptions(m, C, r1, r2, di, only)
% The five circular-orbit models; only: one model's key, or '' for all.
tofmax = m.tof_max_h*3600;
same = abs(r2 - r1) < 1e-6;
want = @(k) isempty(only) || strcmp(only, k);
c = {};

% --- 1. Hohmann, whole plane change at the larger-radius burn -------------
if want('hohmann')
    [dv, tof, nb, plan] = odst_hohmann(r1, r2, di);
    if same, nm = 'Simple plane change'; else, nm = 'Hohmann (combined plane change)'; end
    c{end+1} = newOpt('hohmann', nm, dv, tof, nb, m.isp_chem, C.CPLX.hohmann, plan, struct());
end

% --- 2. Hohmann with numerically optimised plane-change split -------------
if want('hohmann_split')
    s = odst_hohmann_split_opt(r1, r2, di);
    [dv, tof, nb, plan] = odst_hohmann(r1, r2, di, s);
    o = newOpt('hohmann_split', 'Hohmann (optimal plane-change split)', dv, tof, nb, ...
               m.isp_chem, C.CPLX.hohmann_split, plan, struct('split_at_r1', s));
    if di <= 0 || same
        o.applicable = false;
        o.note = 'identical to Hohmann (no plane change to split)';
    end
    c{end+1} = o;
end

% --- 3. Bi-elliptic with optimised intermediate apoapsis rB ---------------
if want('bielliptic')
    rb = odst_bielliptic_opt(r1, r2, di, max(m.rb_max, 1.01*max(r1, r2)), tofmax);
    tofOk = ~isempty(rb);
    if ~tofOk
        rb = 2*max(r1, r2);   % even the fastest bi-elliptic misses the deadline
    end
    [dv, tof, nb, plan] = odst_bielliptic(r1, r2, rb, di);
    o = newOpt('bielliptic', 'Bi-elliptic', dv, tof, nb, m.isp_chem, C.CPLX.bielliptic, ...
               plan, struct('rb', rb, 'rb_over_r1', rb/r1));
    if tofOk && rb < 1.001*max(r1, r2)
        o.applicable = false;
        o.note = 'optimum degenerates to Hohmann (rB -> r2)';
    end
    c{end+1} = o;
end

% --- 4. One-tangent burn (fast transfer) ----------------------------------
if ~same && want('one_tangent')
    [~, tofH] = odst_hohmann(r1, r2, di);
    if tofmax < tofH
        k = odst_onetangent_for_tof(r1, r2, di, tofmax);
        if isempty(k)
            k = 50; note = 'cannot meet TOF limit (k > 50)';
        else
            note = 'k solved to meet TOF limit';
        end
    else
        k = m.k_one_tangent; note = 'default k';
    end
    [dv, tof, nb, plan, info] = odst_onetangent(r1, r2, k, di);
    o = newOpt('one_tangent', 'One-tangent burn (fast)', dv, tof, nb, m.isp_chem, ...
               C.CPLX.one_tangent, plan, info);
    o.note = note;
    c{end+1} = o;
end

% --- 5. Low-thrust spiral (Edelbaum) --------------------------------------
if ~want('low_thrust'), return; end
isp = []; if ~isnan(m.mass0_kg), isp = m.isp_ep; end       % spacecraft: constant thrust, falling mass
[dv, tof] = odst_lowthrust(r1, r2, di, m.accel_ep, isp);
o = newOpt('low_thrust', 'Low-thrust spiral (Edelbaum)', dv, tof, 0, m.isp_ep, ...
           C.CPLX.low_thrust, zeros(0, 3), struct());
if di > 2          % rad, 114.6 deg: Edelbaum's yaw angle sweeps pi/2*di, which must not exceed pi
    o.applicable = false;
    o.note = 'Edelbaum''s closed form holds only up to di = 114.6 deg';
end
c{end+1} = o;
end

% A circular model's burn plan [t, speed after, plane rotation] as a vector plan in the
% relative frame: the burns of odst_simulate, with the coasts solved analytically
% (odst_kepler) instead of integrated, so the plan is exact. Plane changes turn toward the
% target plane, normal [0; -sin(di); cos(di)], and must happen on the node line (x): if the
% first one wouldn't (a raising one-tangent transfer changes plane on arrival), the start
% point is moved back so that it does. Same as circular_vec_plan in odst.py.
function [y0, vplan] = circularVecPlan(mu, r1, plan, di)
y = [r1; 0; 0; 0; sqrt(mu/r1); 0]; t0 = 0; phi = 0;
for i = 1:size(plan, 1)                    % in the start plane up to the first plane change
    if plan(i, 1) > t0, y = odst_kepler(y, plan(i, 1) - t0); t0 = plan(i, 1); end
    if plan(i, 3) ~= 0
        phi = atan2(y(2), y(1));
        if abs(sin(phi)) < 1e-9, phi = 0; end
        break
    end
    y(4:6) = plan(i, 2)*cross([0; 0; 1], y(1:3)/norm(y(1:3)));
end
c = cos(phi); s = sin(phi); v1 = sqrt(mu/r1);
y = [r1*c; -r1*s; 0; v1*s; v1*c; 0]; y0 = y; t0 = 0;
toward = [0; -sin(di); cos(di)];
vplan = zeros(size(plan, 1), 4);
for i = 1:size(plan, 1)
    if plan(i, 1) > t0, y = odst_kepler(y, plan(i, 1) - t0); end
    r = y(1:3); v = y(4:6); nh = cross(r, v)/norm(cross(r, v));
    nA = rodrigues(nh, r, plan(i, 3)); nB = rodrigues(nh, r, -plan(i, 3));
    if dot(nA, toward) >= dot(nB, toward), n = nA; else, n = nB; end
    vNew = plan(i, 2)*cross(n, r/norm(r));
    vplan(i, :) = [plan(i, 1), (vNew - v).'];
    y = [r; vNew]; t0 = plan(i, 1);
end
end

function v = rodrigues(v, axis, ang)
k = axis/norm(axis);
v = v*cos(ang) + cross(k, v)*sin(ang) + k*dot(k, v)*(1 - cos(ang));
end

% plan: circular-model burn plan (odst_simulate); y0 + vplan: vector burn plan of every
% impulsive option (odst_simulate_vec, the GMAT export); legs: [radius, dv] per burn
function o = newOpt(key, name, dv, tof, burns, isp, cplx, plan, params)
o = struct('key', key, 'name', name, 'dv', dv, 'tof', tof, 'burns', burns, ...
           'isp', isp, 'complexity', cplx, 'applicable', true, 'note', '', ...
           'plan', plan, 'y0', [], 'vplan', zeros(0, 4), 'legs', zeros(0, 2), ...
           'params', params, 'dv_loss', 0, 'prop_frac', NaN, 'prop_kg', NaN, 'feasible', true, ...
           'reasons', '', 'chance', NaN, 'pareto', false, 'j2_dv', NaN, 'm_eff', [], 'score_ws', NaN, 'score_topsis', NaN);
end

% ---------------------------------------------------------- J2 (Earth's oblateness)
% With orbital elements and m.j2, J2's secular drift turns every orbit's plane about Earth's axis
% at nodeRate during the transfer. Impulsive transfers are aimed at the target as the spacecraft
% will find it on arrival (j2Plan); the spiral flies to the target's plane as it is, then makes up
% what J2 turned (spiralJ2); nodalDrift uses the drift itself to turn the plane. As in odst.py.
function r = nodeRate(C, a, e, i)
% secular J2 drift of the ascending node [rad/s]: -3/2 n J2 (R/p)^2 cos i
r = -1.5*sqrt(C.MU/a^3)*C.J2*(C.RE/(a*(1 - e*e)))^2*cos(i);
end

function [W1, W2] = normals(m)
e = m.elements*pi/180;
[~, ~, W1] = odst_perifocal(e(1), e(2), e(3));
[~, ~, W2] = odst_perifocal(e(4), e(5), e(6));
end

function [A, B] = tilt(Wa, Wb)
% cos i along the plane turned from normal Wa toward Wb by phi: cos i = A cos(phi) + B sin(phi)
k = cross(Wa, Wb); nk = norm(k);
A = Wa(3); B = 0;
if nk >= 1e-12, B = (k(1)*Wa(2) - k(2)*Wa(1))/nk; end
end

function d = spiralDrift(m, C, ra, rb, Wa, Wb, dv, u0)
% The spacecraft's J2 node drift [rad] along an Edelbaum spiral from radius ra and plane Wa to rb
% and Wb (Simpson in the velocity increment u; dt/du = 1/a, or e^(-u/ve)/a0 at constant thrust,
% after the u0 already flown)
V0 = sqrt(C.MU/ra); V1 = sqrt(C.MU/rb);
di = acos(max(-1, min(1, dot(Wa, Wb))));
b0 = atan2(sin(pi/2*di), V0/V1 - cos(pi/2*di));
[A, B] = tilt(Wa, Wb);
a0 = m.accel_ep*1e-3; ve = m.isp_ep*C.G0; sc = ~isnan(m.mass0_kg);
h = dv/400; s = 0;
for j = 0:400
    u = j*h;
    V = sqrt(V0*V0 - 2*V0*u*cos(b0) + u*u);
    phi = 2/pi*(atan2(V0*sin(b0), V0*cos(b0) - u) - b0);
    a = C.MU/V^2;
    g = 1; if sc, g = exp(-(u0 + u)/ve); end
    w = 2; if j == 0 || j == 400, w = 1; elseif mod(j, 2) == 1, w = 4; end
    s = s + (-1.5*sqrt(C.MU/a^3)*C.J2*(C.RE/a)^2*(A*cos(phi) + B*sin(phi))*g/a0)*w;
end
d = h/3*s;
end

function o = spiralJ2(m, C, o)
% The spiral under J2: Edelbaum's spiral to the target's plane as it is at departure, then
% low-thrust plane changes at the target orbit for what J2 turned the planes apart meanwhile,
% until under 1e-5 rad. If a correction leaves more than the one before, the spiral can't keep
% up with the target's plane, and it isn't offered.
[W1, W2] = normals(m); i2 = m.elements(4)*pi/180; r2 = C.RE + m.alt2;
rn = nodeRate(C, r2, 0, i2); a0 = m.accel_ep*1e-3; ve = m.isp_ep*C.G0;
if isnan(m.mass0_kg), clock = @(u) u/a0; else, clock = @(u) ve/a0*(1 - exp(-u/ve)); end
ra = C.RE + m.alt1; Wa = W1; Wb = W2; u = 0; steps = 0; last = Inf;
while true
    dv = odst_lowthrust(ra, r2, acos(max(-1, min(1, dot(Wa, Wb)))), m.accel_ep);
    t = clock(u + dv) - clock(u);
    lag = rn*t - spiralDrift(m, C, ra, r2, Wa, Wb, dv, u);    % the target's node minus the spacecraft's, over this leg
    u = u + dv; steps = steps + 1;
    left = 2*asin(min(1, sin(i2)*abs(sin(lag/2))));
    if left < 1e-5, break; end
    if left >= last || steps == 30
        o.applicable = false; o.note = 'J2 turns the planes apart faster than the spiral can follow';
        o.dv = NaN; o.tof = NaN;
        return
    end
    ra = r2; Wa = W2; Wb = rodrigues(W2, [0; 0; 1], lag); last = left;
end
o.j2_dv = u - o.dv; o.dv = u; o.tof = clock(u); o.params = struct('corrections', steps - 1);
end

function d = scDrift(m, C, o)
% The spacecraft's J2 node drift [rad] along impulsive option o, planned for mission m: each
% coast at the rate of the orbit it flies, its plane turned by the plane changes so far
[W1, W2] = normals(m); [A, B] = tilt(W1, W2); d = 0;
if ~isempty(o.plan)                 % circular models: the orbit after each burn
    phi = 0;
    for k = 1:size(o.plan, 1) - 1
        t = o.plan(k, 1); v = o.plan(k, 2); r = o.plan(k, 4);
        phi = phi + o.plan(k, 3);
        d = d + nodeRate(C, 1/(2/r - v*v/C.MU), abs(r*v*v/C.MU - 1), ...
                         acos(max(-1, min(1, A*cos(phi) + B*sin(phi)))))*(o.plan(k + 1, 1) - t);
    end
    return
end
y = o.y0(:); t0 = 0; rate = 0;      % elliptical models: the vector plan, turned by m.R
for k = 1:size(o.vplan, 1)
    t = o.vplan(k, 1);
    if t > t0, d = d + rate*(t - t0); y = odst_kepler(y, t - t0); end
    y(4:6) = y(4:6) + o.vplan(k, 2:4).';
    r = y(1:3); v = y(4:6); h = cross(r, v);
    a = 1/(2/norm(r) - dot(v, v)/C.MU); e = norm(cross(v, h)/C.MU - r/norm(r));
    rate = nodeRate(C, a, e, acos(max(-1, min(1, (m.R(3, :)*h)/norm(h))))); t0 = t;
end
end

function me = aimed(m, dn, da)
% The mission with its target turned as the spacecraft will find it on arrival: node by dn,
% periapsis by da [rad]
me = m;
me.elements(5) = m.elements(5) + dn*(180/pi);
me.elements(6) = m.elements(6) + da*(180/pi);
[me.di_deg, me.w1_deg, me.dw_deg, me.R] = odst_relative_geometry(me);
end

function o = buildOne(m, C, key)
% The option `key` for mission m, before the rules
r1 = C.RE + m.alt1; r2 = C.RE + m.alt2; di = m.di_deg*pi/180;
ra1 = r1; ra2 = r2;
if ~isnan(m.apo1), ra1 = C.RE + m.apo1; end
if ~isnan(m.apo2), ra2 = C.RE + m.apo2; end
if ra1 - r1 > 1e-6 || ra2 - r2 > 1e-6
    c = ellipticalOptions(m, C, r1, ra1, r2, ra2, di, key);
else
    c = circularOptions(m, C, r1, r2, di, key);
end
o = c{1};
end

function p = j2Plan(m, C, o)
% Option o planned for J2. An impulsive transfer is aimed at the target as the spacecraft will
% find it on arrival: the node turned by (the target's drift) - (the spacecraft's drift), the
% periapsis by its own drift, iterated to a fixed point (the flight time hardly depends on the
% aim, so it settles; Lambert is searched again once). The spiral: see spiralJ2.
p = o;
if ~o.applicable, return; end
if strcmp(o.key, 'low_thrust'), p = spiralJ2(m, C, o); return; end
r2 = C.RE + m.alt2; ra2 = r2; if ~isnan(m.apo2), ra2 = C.RE + m.apo2; end
ra1 = C.RE + m.alt1; if ~isnan(m.apo1), ra1 = C.RE + m.apo1; end
ell = ra1 - (C.RE + m.alt1) > 1e-6 || ra2 - r2 > 1e-6;
a2 = (r2 + ra2)/2; e2 = (ra2 - r2)/(ra2 + r2); i2 = m.elements(4)*pi/180;
rn = nodeRate(C, a2, e2, i2);
ra = 0.75*sqrt(C.MU/a2^3)*C.J2*(C.RE/(a2*(1 - e2*e2)))^2*(5*cos(i2)^2 - 1);    % the target's apse line
d = rn*o.tof - scDrift(m, C, o); T = o.tof;
n = 50; if ell, n = 1; end
for it = 1:n
    me = aimed(m, d, ra*T);
    p = buildOne(me, C, o.key);
    if ~p.applicable, break; end
    dNew = rn*p.tof - scDrift(me, C, p); Tn = p.tof;
    done = abs(dNew - d) < 1e-10 && abs(Tn - T) <= 1e-9*max(T, 1);
    d = dNew; T = Tn;
    if done, break; end
end
p.m_eff = me; p.j2_dv = p.dv - o.dv;
end

function o = nodalDrift(m, C)
% J2-assisted plane change between circular orbits: Hohmann to a parking orbit of radius rd, wait
% while the two planes' different J2 drift closes the gap between their nodes, then Hohmann to the
% target with the inclination change at its larger-radius burn (on the shared node line). rd is
% the cheapest parking orbit that meets the time limit, found on a grid and refined by bisection to
% the limit; with none, the fastest. As nodal_drift in odst.py.
e = m.elements; i1 = e(1)*pi/180; i2 = e(4)*pi/180; gap0 = (e(5) - e(2))*pi/180;
r1 = C.RE + m.alt1; r2 = C.RE + m.alt2; di = abs(i2 - i1); tofmax = m.tof_max_h*3600;
o = newOpt('nodal_drift', 'Nodal drift (J2)', NaN, NaN, 0, m.isp_chem, C.CPLX.nodal_drift, zeros(0, 5), struct());
if min(sin(i1), sin(i2)) < 1e-6
    o.applicable = false; o.note = 'an equatorial orbit has no node to move'; return
end
if abs(sin(gap0/2)) < 1e-9
    o.applicable = false; o.note = 'the nodes already coincide, so Hohmann does the plane change'; return
end
rn = nodeRate(C, r2, 0, i2);
plan = @(rd) nodalPlan(C, rd, r1, r2, i1, i2, di, gap0, rn);
rlo = C.RE + C.NODAL_ALT_MIN; rhi = 4*max([r1, r2, rlo]);
ends = [r1, r2];
grid = unique([rlo*(rhi/rlo).^((0:399)/399), ends(ends >= rlo)]);
res = zeros(numel(grid), 2);
for k = 1:numel(grid), [res(k, 1), res(k, 2)] = plan(grid(k)); end
fits = @(T) T <= tofmax*(1 + 1e-6);
ok = find(fits(res(:, 2))).';
if isempty(ok)                      % the fastest; of those that tie (a polar orbit's node stays still), the cheapest
    Tmin = min(res(:, 2));
    if isinf(Tmin)
        o.applicable = false; o.note = 'the planes drift at the same rate'; return
    end
    tie = find(res(:, 2) <= Tmin*(1 + 1e-9));
    [~, jj] = min(res(tie, 1)); rd = grid(tie(jj));
else
    [~, jj] = min(res(ok, 1)); j = ok(jj); rd = grid(j);
    for q = [j - 1, j + 1]          % a cheaper neighbour past the time limit: bisect to the limit
        if q >= 1 && q <= numel(grid) && ~any(ok == q) && res(q, 1) < res(j, 1)
            a = grid(j); b = grid(q);
            for it = 1:60
                c = (a + b)/2; [~, Tc] = plan(c);
                if fits(Tc), a = c; else, b = c; end
            end
            if plan(a) < plan(rd), rd = a; end
        end
    end
end
[dv, tof, t1, tw, rho, g] = plan(rd);
[~, ~, ~, p1] = odst_hohmann(r1, rd, 0);
[~, ~, ~, p2] = odst_hohmann(rd, r2, di);
p2(:, 1) = (t1 + tw) + p2(:, 1);
bp = [p1; p2];
o.plan = bp(bp(:, 5) > 1e-9, :);
o.dv = dv; o.tof = tof; o.burns = size(o.plan, 1);
o.params = struct('alt_d', rd - C.RE, 'wait_days', tw/86400, 'rate_deg_day', rho*(180/pi)*86400, ...
                  'gap_deg', g*(180/pi), 'di_deg', di*(180/pi));
end

function [dv, tof, t1, tw, rho, g] = nodalPlan(C, rd, r1, r2, i1, i2, di, gap0, rn)
% The nodal drift through parking orbit rd: dv, flight time, the first transfer's time, the wait,
% the drift rate against the target, the node gap closed while waiting
[dv1, t1] = odst_hohmann(r1, rd, 0);
[dv2, t2] = odst_hohmann(rd, r2, di);
leg = @(a, b, i) nodeRate(C, (a + b)/2, abs(b - a)/(a + b), i);
if r2 >= rd, iLeg = i1; else, iLeg = i2; end
gap = gap0 + rn*(t1 + t2) - leg(r1, rd, i1)*t1 - leg(rd, r2, iLeg)*t2;
rho = nodeRate(C, rd, 0, i1) - rn;
g = mod(gap, 2*pi);
if min(g, 2*pi - g) < 1e-12, g = 0; elseif rho <= 0, g = g - 2*pi; end
if g == 0, tw = 0; elseif abs(rho) > 1e-18, tw = g/rho; else, tw = Inf; end
dv = dv1 + dv2; tof = t1 + tw + t2;
end

% ------------------------------------------- rule-based feasibility filter
function opts = applyRules(m, opts)
tofmax = m.tof_max_h*3600;
for i = 1:numel(opts)
    o = opts(i);
    r = {};
    o.feasible = o.applicable;
    if ~o.applicable
        r{end+1} = o.note; %#ok<AGROW>
    else
        need = (o.dv + o.dv_loss)*(1 + m.margin_pct/100);
        if need > m.dv_budget
            o.feasible = false;
            r{end+1} = sprintf('dv %.3f + loss %.3f + %g%% margin = %.3f > budget %.3f km/s', ...
                               o.dv, o.dv_loss, m.margin_pct, need, m.dv_budget); %#ok<AGROW>
        end
        if o.tof > tofmax*(1 + 1e-6)
            o.feasible = false;
            r{end+1} = sprintf('TOF %.1f h > limit %.1f h', o.tof/3600, m.tof_max_h); %#ok<AGROW>
        end
        if strcmp(o.key, 'low_thrust') && strcmp(m.propulsion, 'chemical')
            o.feasible = false;
            r{end+1} = 'no electric propulsion on board'; %#ok<AGROW>
        end
        if ~strcmp(o.key, 'low_thrust') && strcmp(m.propulsion, 'electric')
            o.feasible = false;
            r{end+1} = 'no high-thrust (chemical) propulsion on board'; %#ok<AGROW>
        end
        if ~isnan(m.mass0_kg)          % the propellant on board must cover the same dv as the budget
            C = odst_const();
            if strcmp(o.key, 'low_thrust'), onboard = m.prop_ep_kg; else, onboard = m.prop_chem_kg; end
            prop = m.mass0_kg*(1 - exp(-need/(o.isp*C.G0)));
            if prop > onboard
                o.feasible = false;
                r{end+1} = sprintf('needs %.0f kg of propellant with the %g%% margin, %.0f kg on board', ...
                                   prop, m.margin_pct, onboard); %#ok<AGROW>
            end
        end
        if m.uncertainty && o.feasible     % out if it fits the limits in under the required share of the runs
            o.chance = chance(m, o);
            if o.chance < m.confidence_pct/100
                o.feasible = false;
                r{end+1} = sprintf('fits the limits in %d%% of runs with the errors, %g%% required', ...
                                   floor(100*o.chance + 1e-9), m.confidence_pct); %#ok<AGROW>
            end
        end
    end
    o.reasons = strjoin(r, '; ');
    opts(i) = o;
end
end

% Share of the Monte Carlo runs (odst_mc_normals) in which option o fits the budget, the time
% limit and the propellant on board (chance in odst.py). Each run draws an Isp and a thrust
% error and, for each impulsive burn, a pointing error, trimmed straight after the burn at a
% cost of |dv flown - dv planned| = 2 dv sin(theta/2). Burns end at their planned dv
% (accelerometer cutoff), so the thrust error changes the gravity loss and the low-thrust time.
function p = chance(m, o)
C = odst_const(); Z = odst_mc_normals();
isp = o.isp*(1 + m.sigma_isp_pct/100*Z(:, 1));
f = 1 + m.sigma_thrust_pct/100*Z(:, 2);
sc = ~isnan(m.mass0_kg);
if strcmp(o.key, 'low_thrust')
    dv = o.dv; a = m.accel_ep*1e-3*f; ve = isp*C.G0;
    if sc, tof = ve./a.*(1 - exp(-dv./ve)); else, tof = dv./a; end
else
    th = m.sigma_point_deg*pi/180; trims = 0;
    for j = 1:size(o.legs, 1)
        trims = trims + 2*o.legs(j, 2)*sin(th*hypot(Z(:, 1 + 2*j), Z(:, 2 + 2*j))/2);
    end
    dv = o.dv + odst_gravity_loss(o.legs, isp, m.tw_chem*f) + trims;
    tof = o.tof;
end
need = dv*(1 + m.margin_pct/100);
ok = need <= m.dv_budget & tof <= m.tof_max_h*3600*(1 + 1e-6);
if sc
    if strcmp(o.key, 'low_thrust'), onboard = m.prop_ep_kg; else, onboard = m.prop_chem_kg; end
    ok = ok & m.mass0_kg*(1 - exp(-need./(isp*C.G0))) <= onboard;
end
p = mean(ok);
end

% The Pareto front: feasible options that no feasible option beats on both propellant and time
function opts = paretoFront(opts)
f = find([opts.feasible]);
for i = f
    o = opts(i); beaten = false;
    for j = f
        p = opts(j);
        beaten = beaten || (p.prop_frac <= o.prop_frac*(1 + 1e-9) && p.tof <= o.tof*(1 + 1e-9) && ...
                            (p.prop_frac < o.prop_frac*(1 - 1e-9) || p.tof < o.tof*(1 - 1e-9)));
    end
    opts(i).pareto = ~beaten;
end
end

% ------------------------------------------------ multi-criteria scoring
function opts = scoreOptions(m, opts)
idx = find([opts.feasible]);
if isempty(idx), return; end
w = m.weights(:).'/sum(m.weights);
X = zeros(numel(idx), 4);
for j = 1:numel(idx)
    o = opts(idx(j));
    % TOF on a log scale: hours vs months differ in kind, not linearly
    X(j, :) = [o.dv, o.prop_frac, log10(max(o.tof, 60)), o.complexity];
end
% (1) weighted sum, min-max normalised; all criteria are costs -> 1 = best
lo = min(X, [], 1); hi = max(X, [], 1); span = hi - lo;
U = ones(size(X));
for k = 1:4
    if span(k) > 1e-12
        U(:, k) = (hi(k) - X(:, k))/span(k);
    end
end
ws = U*w.';
% (2) TOPSIS, vector normalisation
nrm = sqrt(sum(X.^2, 1)); nrm(nrm < 1e-12) = 1;
V = X./nrm.*w;
best = min(V, [], 1); worst = max(V, [], 1);
dp = sqrt(sum((V - best).^2, 2));
dm = sqrt(sum((V - worst).^2, 2));
cc = ones(size(dp));
ok = (dp + dm) > 1e-12;
cc(ok) = dm(ok)./(dp(ok) + dm(ok));
for j = 1:numel(idx)
    opts(idx(j)).score_ws = ws(j);
    opts(idx(j)).score_topsis = cc(j);
end
end

% ------------------------------------------ rule-based decision matrix
function adv = ruleAdvice(m, r1, r2, di)
R = max(r1, r2)/min(r1, r2);
adv = {};
if R < 11.94
    adv{end+1} = sprintf('R = %.2f < 11.94: without a plane change, Hohmann beats every bi-elliptic transfer.', R);
elseif R < 15.58
    adv{end+1} = sprintf('11.94 <= R = %.2f < 15.58: without a plane change, bi-elliptic wins only if rB is large enough.', R);
else
    adv{end+1} = sprintf('R = %.2f >= 15.58: without a plane change, bi-elliptic beats Hohmann for any rB > r2.', R);
end
if m.di_deg > 0
    adv{end+1} = ['Plane change: perform it where speed is lowest (largest radius) or combine with a burn. ' ...
                  'That can make bi-elliptic the cheapest at any R.'];
end
if R < 1 + 1e-9 && m.di_deg > 48.94
    adv{end+1} = 'di > 48.94 deg: a three-impulse (bi-elliptic) plane change beats a single burn.';
end
[~, tofH] = odst_hohmann(r1, r2, di);
if m.tof_max_h < tofH/3600
    adv{end+1} = sprintf(['TOF limit %.2f h < Hohmann %.2f h: a fast (one-tangent) ' ...
                          'transfer is required.'], m.tof_max_h, tofH/3600);
end
if ~strcmp(m.propulsion, 'chemical') && di <= 2       % Edelbaum's limit, 114.6 deg
    adv{end+1} = 'Electric propulsion available: low-thrust spiral minimises propellant if time allows.';
end
end
