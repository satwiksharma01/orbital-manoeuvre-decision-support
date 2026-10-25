function txt = odst_burns_json(m, o)
%ODST_BURNS_JSON Option o's burn plan as JSON: the start state, then each burn's time and
%   inertial dv (km, s, km/s). Same fields as burns_json in odst.py.
%   txt = odst_burns_json(m, o)          m: odst_mission, o: an option from odst_evaluate
C = odst_const();
a1 = m.alt1; if ~isnan(m.apo1), a1 = m.apo1; end
a2 = m.alt2; if ~isnan(m.apo2), a2 = m.apo2; end
b = cell(1, size(o.vplan, 1));
for k = 1:numel(b)
    dv = o.vplan(k, 2:4) + 0;            % + 0 turns -0 into 0
    b{k} = struct('t_s', o.vplan(k, 1), 'dv_kms', dv, 'dv_mag_kms', norm(dv));
end
s = struct('tool', 'Orbital Manoeuvre Decision Support Tool', 'transfer', o.name, 'dv_kms', o.dv, ...
           'tof_s', o.tof, 'start_orbit_alt_km', [m.alt1 a1], 'target_orbit_alt_km', [m.alt2 a2], ...
           'di_deg', m.di_deg, 'dw_deg', m.dw_deg, 'frame', C.FRAME, 'epoch_utc', C.EPOCH, ...
           'mu_km3_s2', C.MU, 'start_state_km_kms', o.y0(:).' + 0, 'burns', {b});
txt = [jsonencode(s, 'PrettyPrint', true) newline];
end
