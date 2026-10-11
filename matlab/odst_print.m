function odst_print(m, opts, rec, advice)
%ODST_PRINT Print the decision table, recommendation and rule advice.
fprintf('\nMission: alt %.0f -> %.0f km, di = %.2f deg, budget %.2f km/s, TOF max %.1f h\n', ...
        m.alt1, m.alt2, m.di_deg, m.dv_budget, m.tof_max_h);
fprintf('Weights [dv prop TOF cplx] = [%s], propulsion: %s\n', num2str(m.weights), m.propulsion);
fprintf('%-38s %8s %11s %5s %7s %5s %6s %7s  %s\n', 'Option', 'dv km/s', 'TOF', ...
        'burns', 'mp/m0', 'cplx', 'WS', 'TOPSIS', 'status');
fprintf('%s\n', repmat('-', 1, 104));
for i = 1:numel(opts)
    o = opts(i);
    if o.feasible
        st = 'OK';
        ws = sprintf('%.3f', o.score_ws); tp = sprintf('%.3f', o.score_topsis);
    else
        st = ['X ' o.reasons]; ws = '-'; tp = '-';
    end
    fprintf('%-38s %8.4f %11s %5d %7.3f %5.1f %6s %7s  %s\n', o.name, o.dv, ...
            odst_fmt_time(o.tof), o.burns, o.prop_frac, o.complexity, ws, tp, st);
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
