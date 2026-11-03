function odst_print(m, opts, rec, advice)
%ODST_PRINT Print the decision table, recommendation and rule advice.
o1 = sprintf('%.0f', m.alt1); o2 = sprintf('%.0f', m.alt2); dw = '';
if ~isnan(m.apo1), o1 = sprintf('%.0f x %.0f', m.alt1, m.apo1); end
if ~isnan(m.apo2), o2 = sprintf('%.0f x %.0f', m.alt2, m.apo2); end
if ~isnan(m.apo1) || ~isnan(m.apo2), dw = sprintf(', dw = %g deg', m.dw_deg); end
fprintf('\nMission: alt %s -> %s km, di = %.2f deg%s, budget %.2f km/s (margin %g%%, T/W %g), TOF max %.1f h\n', ...
        o1, o2, m.di_deg, dw, m.dv_budget, m.margin_pct, m.tw_chem, m.tof_max_h);
fprintf('Weights [dv prop TOF cplx] = [%s], propulsion: %s\n', num2str(m.weights), m.propulsion);
fit = '';
if m.uncertainty
    fprintf('Uncertainty, 1 sigma: Isp %g%%, thrust %g%%, pointing %g deg per axis; 4000 runs, %g%% required\n', ...
            m.sigma_isp_pct, m.sigma_thrust_pct, m.sigma_point_deg, m.confidence_pct);
    fit = sprintf(' %5s', 'fit');
end
j2 = m.j2 && ~isempty(m.elements); j2c = '';
if j2
    fprintf('J2: each plan allows for the planes'' drift; ''J2 m/s'' is its dv minus the two-body plan''s\n');
    j2c = sprintf(' %7s', 'J2 m/s');
end
fprintf('%-38s %8s %8s %11s %5s %7s %5s%s%s %6s %7s  %s\n', 'Option', 'dv km/s', 'loss m/s', 'TOF', ...
        'burns', 'mp/m0', 'cplx', fit, j2c, 'WS', 'TOPSIS', 'status');
fprintf('%s\n', repmat('-', 1, 113 + numel(fit) + numel(j2c)));
for i = 1:numel(opts)
    o = opts(i);
    if m.uncertainty
        if isnan(o.chance), fit = sprintf(' %5s', '-'); else, fit = sprintf(' %4d%%', floor(100*o.chance + 1e-9)); end
    end
    if o.feasible
        st = 'OK';
        if o.pareto, st = 'OK, Pareto front'; end
        ws = sprintf('%.3f', o.score_ws); tp = sprintf('%.3f', o.score_topsis);
    else
        st = ['X ' o.reasons]; ws = '-'; tp = '-';
    end
    if j2
        if isnan(o.j2_dv), j2c = sprintf(' %7s', '-'); else, j2c = sprintf(' %+7.1f', 1e3*o.j2_dv); end
    end
    fprintf('%-38s %8.4f %8.1f %11s %5d %7.3f %5.1f%s%s %6s %7s  %s\n', o.name, o.dv, 1e3*o.dv_loss, ...
            odst_fmt_time(o.tof), o.burns, o.prop_frac, o.complexity, fit, j2c, ws, tp, st);
    if strcmp(o.key, 'nodal_drift') && o.applicable
        fprintf('%-38s wait %.2f days at %.0f km, where the planes'' nodes close at %+.3f deg/day (%+.2f deg)\n', '', ...
                o.params.wait_days, o.params.alt_d, o.params.rate_deg_day, o.params.gap_deg);
    end
end
if isempty(rec)
    fprintf('Recommended: NONE FEASIBLE - relax the constraints\n');
else
    fprintf('Recommended: %s\n', rec.name);
end
for i = 1:numel(advice)
    fprintf('  rule: %s\n', advice{i});
end
end
