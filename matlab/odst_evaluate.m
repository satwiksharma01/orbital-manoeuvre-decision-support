function [opts, rec, advice] = odst_evaluate(m)
%ODST_EVALUATE Evaluate all candidate manoeuvres for mission m and rank them.
%   [opts, rec, advice] = odst_evaluate(odst_mission(...))
%
%   opts   struct array, one element per candidate manoeuvre:
%          key, name, dv [km/s], tof [s], burns, isp [s], complexity,
%          applicable, note, plan (burn plan), params, prop_frac (mp/m0),
%          feasible, reasons, score_ws, score_topsis
%   rec    recommended option (highest weighted-sum score among feasible), or []
%   advice cellstr of rule-of-thumb statements (rule-based decision matrix)
%
%   Pipeline:  1) analytical models  2) numerical optimisation of free
%   parameters  3) rule-based feasibility filter  4) multi-criteria scoring.
%   If either orbit is elliptical (apo1/apo2 set), the apse-to-apse and
%   Lambert transfers are sized instead of the five circular-orbit models.
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
    c = ellipticalOptions(m, C, r1, ra1, r2, ra2, di);
else
    c = circularOptions(m, C, r1, r2, di);
end

opts = [c{:}];
for i = 1:numel(opts)                                   % Tsiolkovsky
    opts(i).prop_frac = 1 - exp(-opts(i).dv/(opts(i).isp*C.G0));
    if ~isempty(opts(i).plan)
        opts(i).legs = opts(i).plan(:, 4:5);
    end
    if ~isempty(opts(i).legs)
        opts(i).dv_loss = odst_gravity_loss(opts(i).legs, opts(i).isp, m.tw_chem);
    end
end
opts = applyRules(m, opts);
opts = scoreOptions(m, opts);

rec = [];
f = find([opts.feasible]);
if ~isempty(f)
    [~, j] = max([opts(f).score_ws]);
    rec = opts(f(j));
end
if ell
    advice = ellipticalAdvice(m, opts(1).applicable);
else
    advice = ruleAdvice(m, r1, r2, di);
end
end

% =========================================================================
function c = ellipticalOptions(m, C, r1, ra1, r2, ra2, di)
% Hohmann, bi-elliptic, one-tangent and Edelbaum assume circular orbits, so an
% elliptical mission is sized with the apse-to-apse and Lambert models only.
dw = m.dw_deg*pi/180;
c = {};
[res, why] = odst_apse(r1, ra1, r2, ra2, di, dw);
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
res = odst_lambert_opt(r1, ra1, r2, ra2, di, dw, m.tof_max_h*3600);
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

function adv = ellipticalAdvice(m, apseOk)
adv = {['Elliptical orbits: Hohmann, bi-elliptic, one-tangent and Edelbaum low thrust assume ' ...
        'circular orbits, so the apse-to-apse and Lambert transfers are sized instead.']};
if apseOk
    adv{end+1} = 'Apse lines aligned: a transfer between apsides is available, usually the cheapest.';
else
    adv{end+1} = sprintf('Apse lines %g deg apart: only a Lambert transfer is offered.', m.dw_deg);
end
end

function c = circularOptions(m, C, r1, r2, di)
tofmax = m.tof_max_h*3600;
same = abs(r2 - r1) < 1e-6;
c = {};

% --- 1. Hohmann, whole plane change at the larger-radius burn -------------
[dv, tof, nb, plan] = odst_hohmann(r1, r2, di);
if same, nm = 'Simple plane change'; else, nm = 'Hohmann (combined plane change)'; end
c{end+1} = newOpt('hohmann', nm, dv, tof, nb, m.isp_chem, C.CPLX.hohmann, plan, struct());

% --- 2. Hohmann with numerically optimised plane-change split -------------
s = odst_hohmann_split_opt(r1, r2, di);
[dv, tof, nb, plan] = odst_hohmann(r1, r2, di, s);
o = newOpt('hohmann_split', 'Hohmann (optimal plane-change split)', dv, tof, nb, ...
           m.isp_chem, C.CPLX.hohmann_split, plan, struct('split_at_r1', s));
if di <= 0 || same
    o.applicable = false;
    o.note = 'identical to Hohmann (no plane change to split)';
end
c{end+1} = o;

% --- 3. Bi-elliptic with optimised intermediate apoapsis rB ---------------
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

% --- 4. One-tangent burn (fast transfer) ----------------------------------
if ~same
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
[dv, tof] = odst_lowthrust(r1, r2, di, m.accel_ep);
c{end+1} = newOpt('low_thrust', 'Low-thrust spiral (Edelbaum)', dv, tof, 0, m.isp_ep, ...
                  C.CPLX.low_thrust, zeros(0, 3), struct());
end

% plan: circular-model burn plan (odst_simulate); y0 + vplan: vector burn plan
% of the elliptical models (odst_simulate_vec); legs: [radius, dv] per burn
function o = newOpt(key, name, dv, tof, burns, isp, cplx, plan, params)
o = struct('key', key, 'name', name, 'dv', dv, 'tof', tof, 'burns', burns, ...
           'isp', isp, 'complexity', cplx, 'applicable', true, 'note', '', ...
           'plan', plan, 'y0', [], 'vplan', zeros(0, 4), 'legs', zeros(0, 2), ...
           'params', params, 'dv_loss', 0, 'prop_frac', NaN, 'feasible', true, ...
           'reasons', '', 'score_ws', NaN, 'score_topsis', NaN);
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
    end
    o.reasons = strjoin(r, '; ');
    opts(i) = o;
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
    adv{end+1} = sprintf('R = %.2f < 11.94: Hohmann is the minimum-dv two/three-impulse transfer.', R);
elseif R < 15.58
    adv{end+1} = sprintf('11.94 <= R = %.2f < 15.58: bi-elliptic wins only if rB is large enough.', R);
else
    adv{end+1} = sprintf('R = %.2f >= 15.58: bi-elliptic beats Hohmann for any rB > r2.', R);
end
if m.di_deg > 0
    adv{end+1} = 'Plane change: perform it where speed is lowest (largest radius) or combine with a burn.';
end
if R < 1 + 1e-9 && m.di_deg > 48.94
    adv{end+1} = 'di > 48.94 deg: a three-impulse (bi-elliptic) plane change beats a single burn.';
end
[~, tofH] = odst_hohmann(r1, r2, di);
if m.tof_max_h < tofH/3600
    adv{end+1} = sprintf(['TOF limit %.2f h < Hohmann %.2f h: a fast (one-tangent) ' ...
                          'transfer is required.'], m.tof_max_h, tofH/3600);
end
if ~strcmp(m.propulsion, 'chemical')
    adv{end+1} = 'Electric propulsion available: low-thrust spiral minimises propellant if time allows.';
end
end
