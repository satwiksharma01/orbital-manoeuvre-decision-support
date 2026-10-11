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
C = odst_const();
r1 = C.RE + m.alt1;
r2 = C.RE + m.alt2;
di = m.di_deg*pi/180;
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

opts = [c{:}];
for i = 1:numel(opts)                                   % Tsiolkovsky
    opts(i).prop_frac = 1 - exp(-opts(i).dv/(opts(i).isp*C.G0));
end
opts = applyRules(m, opts);
opts = scoreOptions(m, opts);

rec = [];
f = find([opts.feasible]);
if ~isempty(f)
    [~, j] = max([opts(f).score_ws]);
    rec = opts(f(j));
end
advice = ruleAdvice(m, r1, r2, di);
end

% =========================================================================
function o = newOpt(key, name, dv, tof, burns, isp, cplx, plan, params)
o = struct('key', key, 'name', name, 'dv', dv, 'tof', tof, 'burns', burns, ...
           'isp', isp, 'complexity', cplx, 'applicable', true, 'note', '', ...
           'plan', plan, 'params', params, 'prop_frac', NaN, 'feasible', true, ...
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
        if o.dv > m.dv_budget
            o.feasible = false;
            r{end+1} = sprintf('dv %.3f > budget %.3f km/s', o.dv, m.dv_budget); %#ok<AGROW>
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
