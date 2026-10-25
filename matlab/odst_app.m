function odst_app()
%ODST_APP Orbital Manoeuvre Planner - interactive decision-support GUI.
%   Same design as the web planner (web/index.html) and the report: example
%   missions in the top bar, grouped mission inputs on the left, and the
%   recommendation, trade space, transfer geometry, ranking and rule-based
%   checks on the right. Every input change re-evaluates all candidate
%   manoeuvres; picking an example or a table row draws the transfer with a
%   critically damped spring timing, like the web version.
%   Requires MATLAB R2021a or newer (uifigure, uigridlayout, uistyle).

T = odst_theme();
C = odst_const();
[pnames, W] = odst_presets();
ARW = char(8594); DEG = char(176); DEL = char(916);
CHK = char(10003); WRN = char(9888); DOT = char(9679); BUL = char(8226); MID = char(183);

defaults = struct('alt1', 400, 'alt2', 35786, 'apo1', NaN, 'apo2', NaN, 'dw_deg', 0, 'di_deg', 28.5, 'dv_budget', 6.5, 'margin_pct', 5, 'tof_max_h', 4800, ...
                  'rb_max', 384400, 'propulsion', 'both', 'isp_chem', 320, 'tw_chem', 0.5, 'isp_ep', 1600, ...
                  'accel_ep', 5e-4, 'preset', 'Balanced', 'elements', []);
defEl = [28.5 0 0 0 0 0];       % element inputs before any are entered: LEO at 28.5 deg -> GEO
examples = {
    ['LEO ' ARW ' GEO'],         struct()
    'High orbit, R = 15',       struct('alt1', 621.863, 'alt2', 98621.863, 'di_deg', 0, 'tof_max_h', 1440, ...
                                       'propulsion', 'chemical', 'preset', 'Propellant-critical')
    '36-min rapid response',    struct('alt1', 400, 'alt2', 800, 'di_deg', 0, 'dv_budget', 1, ...
                                       'tof_max_h', 0.6, 'preset', 'Time-critical')
    'All-electric GEO',         struct('tof_max_h', 8760, 'preset', 'Propellant-critical')
    ['60' DEG ' plane change'], struct('alt1', 500, 'alt2', 500, 'di_deg', 60, 'dv_budget', 10, ...
                                       'tof_max_h', 4800, 'propulsion', 'chemical', 'preset', 'Propellant-critical')
    ['GTO ' ARW ' GEO'],         struct('alt1', 250, 'apo1', 35786, 'alt2', 35786, 'propulsion', 'chemical')
    ['LEO ' ARW ' Molniya'],     struct('alt2', 600, 'apo2', 39700, 'dw_deg', 270, 'di_deg', 0, ...
                                       'propulsion', 'chemical')
    'Kourou GTO',               struct('alt1', 250, 'apo1', 35786, 'elements', [6 0 178 0 0 0], ...
                                       'propulsion', 'chemical')};

% ---- state shared by the nested callbacks
prop = 'both'; shape1 = 'circular'; shape2 = 'circular'; orient = 'relative'; selected = ''; manual = false; order = {};
lastM = []; lastOpts = []; lastRec = []; lastAdv = {};
h = struct();

% ================================================================ window
fig = uifigure('Name', 'Orbital Manoeuvre Planner', 'Position', [40 40 1380 900]);
if isprop(fig, 'Theme'), fig.Theme = 'light'; end
fig.Color = T.bg;
root = uigridlayout(fig, [2 1], 'RowHeight', {58, '1x'}, 'Padding', [0 0 0 0], ...
                    'RowSpacing', 0, 'BackgroundColor', T.bg);

% ---- top bar: brand + example missions
nEx = size(examples, 1);
topBar = uigridlayout(root, [1 4 + nEx], 'Padding', [20 12 20 12], 'ColumnSpacing', 8, ...
                      'BackgroundColor', T.surface);
topBar.ColumnWidth = [{'fit', 'fit', '1x', 'fit'}, repmat({'fit'}, 1, nEx)];
lbl(topBar, 'Orbital Manoeuvre Planner', 16, T.label, 'bold');
lbl(topBar, ['ASE3004 ' MID ' Project 38 ' MID ' Satwik Sharma'], 12, T.label2);
lbl(topBar, '', 12, T.label2);
lbl(topBar, 'Try a mission', 12.5, T.label2);
h.chips = gobjects(1, nEx);
for i = 1:nEx
    h.chips(i) = uibutton(topBar, 'Text', examples{i, 1}, 'FontName', T.font, 'FontSize', 12.5, ...
        'BackgroundColor', T.fill, 'FontColor', T.label, 'ButtonPushedFcn', @(~, ~) loadExample(i));
end

body = uigridlayout(root, [1 2], 'ColumnWidth', {340, '1x'}, 'Padding', [20 16 20 16], ...
                    'ColumnSpacing', 20, 'BackgroundColor', T.bg);

% ================================================================ inputs (left)
left = uigridlayout(body, [10 1], 'Padding', [0 0 0 0], 'RowSpacing', 6, 'BackgroundColor', T.bg);
left.RowHeight = {18, 128, 22, 128, 22, 38, 22, 208, 22, 208};
if isprop(left, 'Scrollable'), left.Scrollable = 'on'; end

groupTitle(left, 'Orbits');
h.orbGrid = listPanel(left, 16);
h.orientSeg = segRow(h.orbGrid, 'Orientation', {[DEL 'i'], 'Elements'}, {'relative', 'elements'}, @pickOrient);
h.shape1 = shapeRow(h.orbGrid, 'Initial orbit', 1);
[h.alt1, h.alt1Lbl] = numRow(h.orbGrid, 'Initial altitude', defaults.alt1, [0 Inf], 'km');
h.apo1 = numRow(h.orbGrid, 'Initial apoapsis', defaults.alt1, [0 Inf], 'km');
h.el = gobjects(1, 6);
h.el(1) = numRow(h.orbGrid, 'Initial inclination', defEl(1), [0 180], 'deg');
h.el(2) = numRow(h.orbGrid, 'Initial RAAN', defEl(2), [0 360], 'deg');
h.el(3) = numRow(h.orbGrid, 'Initial arg. of periapsis', defEl(3), [0 360], 'deg');
h.shape2 = shapeRow(h.orbGrid, 'Target orbit', 2);
[h.alt2, h.alt2Lbl] = numRow(h.orbGrid, 'Target altitude', defaults.alt2, [0 Inf], 'km');
h.apo2 = numRow(h.orbGrid, 'Target apoapsis', defaults.alt2, [0 Inf], 'km');
h.el(4) = numRow(h.orbGrid, 'Target inclination', defEl(4), [0 180], 'deg');
h.el(5) = numRow(h.orbGrid, 'Target RAAN', defEl(5), [0 360], 'deg');
h.el(6) = numRow(h.orbGrid, 'Target arg. of periapsis', defEl(6), [0 360], 'deg');
h.dw   = numRow(h.orbGrid, ['Target apse angle ' DEL char(969)], defaults.dw_deg, [0 360], 'deg');
h.di   = numRow(h.orbGrid, 'Inclination change', defaults.di_deg, [0 180], 'deg');
h.planes = lbl(h.orbGrid, '', 12.5, T.label2);

groupTitle(left, 'Constraints');
g = listPanel(left, 3);
h.dvb = numRow(g, [DEL 'v budget'], defaults.dv_budget, [0 Inf], 'km/s');
h.mrg = numRow(g, [DEL 'v margin'], defaults.margin_pct, [0 100], '%');
h.tof = numRow(g, 'Maximum transfer time', defaults.tof_max_h, [0 Inf], 'h');

groupTitle(left, 'Propulsion on board');
h.propSeg = segmented(left, {'Both', 'Chemical', 'Electric'}, {'both', 'chemical', 'electric'}, @pickProp);

pt = uigridlayout(left, [1 2], 'ColumnWidth', {'1x', 'fit'}, 'Padding', [0 0 0 0], 'BackgroundColor', T.bg);
lbl(pt, '  Mission priorities', 12.5, T.label2, 'bold');
h.customNote = lbl(pt, 'Custom  ', 12.5, T.accent, 'bold');
h.customNote.Visible = 'off';
pp = uipanel(left, 'BorderType', 'none', 'BackgroundColor', T.surface);
pg = uigridlayout(pp, [5 1], 'RowHeight', {36, 40, 40, 40, 40}, 'Padding', [10 8 10 4], ...
                  'RowSpacing', 0, 'BackgroundColor', T.surface);
h.presetSeg = segmented(pg, {'Balanced', 'Propellant', 'Time', 'Simple'}, pnames, @pickPreset);
wNames = {[DEL 'v'], 'Propellant', 'Transfer time', 'Simplicity'};
h.w = gobjects(1, 4); h.wLbl = gobjects(1, 4);
for k = 1:4
    [h.w(k), h.wLbl(k)] = sliderRow(pg, wNames{k}, W(1, k), k);
end

groupTitle(left, 'Engine and limits');
g = listPanel(left, 5);
h.ispc = numRow(g, 'Chemical Isp', defaults.isp_chem, [1 Inf], 's');
h.tw   = numRow(g, 'Chemical thrust/weight', defaults.tw_chem, [0.01 Inf], '');
h.ispe = numRow(g, 'Electric Isp', defaults.isp_ep, [1 Inf], 's');
h.acc  = numRow(g, 'Electric thrust F/m', defaults.accel_ep, [1e-7 1], ['m/s' char(178)]);
h.rbm  = numRow(g, 'Max bi-elliptic apoapsis', defaults.rb_max, [1 Inf], 'km');

% ================================================================ results (right)
right = uigridlayout(body, [3 1], 'RowHeight', {214, '1x', 252}, 'Padding', [0 0 0 0], ...
                     'RowSpacing', 16, 'BackgroundColor', T.bg);

% ---- recommendation card
vp = uipanel(right, 'BorderType', 'none', 'BackgroundColor', T.surface);
vg = uigridlayout(vp, [5 1], 'RowHeight', {22, 38, 22, 54, '1x'}, 'Padding', [24 14 24 12], ...
                  'RowSpacing', 2, 'BackgroundColor', T.surface);
vh = uigridlayout(vg, [1 4], 'ColumnWidth', {16, 'fit', '1x', 'fit'}, 'Padding', [0 0 0 0], ...
                  'ColumnSpacing', 6, 'BackgroundColor', T.surface);
h.vDot = lbl(vh, DOT, 16, T.accent);
h.vEyebrow = lbl(vh, 'Recommended manoeuvre', 13, T.label2, 'bold');
lbl(vh, '', 12, T.label2);
h.vChip = lbl(vh, '', 12, T.good, 'bold');
h.vChip.HorizontalAlignment = 'center';
h.vName = lbl(vg, '', 30, T.label, 'bold');
h.vSub = lbl(vg, '', 15, T.label2);
st = uigridlayout(vg, [2 5], 'RowHeight', {30, 18}, 'ColumnWidth', repmat({'1x'}, 1, 5), ...
                  'Padding', [0 0 0 0], 'RowSpacing', 0, 'ColumnSpacing', 12, 'BackgroundColor', T.surface);
statNames = {['Total ' DEL 'v'], 'Transfer time', 'Burns', 'Propellant (% of mass)', 'Score (of 1)'};
h.stat = gobjects(1, 5);
for k = 1:5
    h.stat(k) = lbl(st, '', 22, T.label, 'bold');
    h.stat(k).Layout.Row = 1; h.stat(k).Layout.Column = k;
    c = lbl(st, statNames{k}, 12, T.label2);
    c.Layout.Row = 2; c.Layout.Column = k;
end
h.vWhy = lbl(vg, '', 13, T.label2);
h.vWhy.WordWrap = 'on'; h.vWhy.VerticalAlignment = 'top';

% ---- trade space + transfer geometry
charts = uigridlayout(right, [1 2], 'ColumnWidth', {'1.25x', '1x'}, 'ColumnSpacing', 16, ...
                      'Padding', [0 0 0 0], 'BackgroundColor', T.bg);
[h.axTrade, ~, ~] = chartCard(charts, 'Trade space', [DEL 'v against transfer time. Hollow markers break a constraint.']);
[h.axGeo, h.geoSub, h.geoTitle, h.export] = chartCard(charts, 'Transfer geometry', '');
h.export.Text = 'Export to GMAT';
h.export.Tooltip = 'Save a NASA GMAT script and a JSON file of the burns';
h.export.ButtonPushedFcn = @(~, ~) exportPlan();

% ---- ranking + rule-based checks
bottom = uigridlayout(right, [1 2], 'ColumnWidth', {'2.3x', '1x'}, 'ColumnSpacing', 16, ...
                      'Padding', [0 0 0 0], 'BackgroundColor', T.bg);
rp = uipanel(bottom, 'BorderType', 'none', 'BackgroundColor', T.surface);
rg = uigridlayout(rp, [3 1], 'RowHeight', {22, 18, '1x'}, 'Padding', [18 12 18 12], ...
                  'RowSpacing', 2, 'BackgroundColor', T.surface);
lbl(rg, 'Ranking', 16, T.label, 'bold');
lbl(rg, 'Select a row to draw that manoeuvre. Scores run from 0 to 1 across the feasible options.', 12, T.label2);
h.tbl = uitable(rg, 'ColumnName', {'Manoeuvre', [DEL 'v km/s'], 'Time', 'Burns', 'Propellant', ...
                'Complexity', 'Score', 'TOPSIS', 'Status'}, 'RowName', {}, 'FontName', T.font, ...
                'FontSize', 12.5, 'CellSelectionCallback', @(~, e) rowPicked(e));
h.tbl.BackgroundColor = T.surface;
if isprop(h.tbl, 'RowStriping'), h.tbl.RowStriping = 'off'; end
if isprop(h.tbl, 'SelectionType'), h.tbl.SelectionType = 'row'; end

cp = uipanel(bottom, 'BorderType', 'none', 'BackgroundColor', T.surface);
cg = uigridlayout(cp, [3 1], 'RowHeight', {22, 18, '1x'}, 'Padding', [18 12 18 12], ...
                  'RowSpacing', 2, 'BackgroundColor', T.surface);
lbl(cg, 'Rule-based checks', 16, T.label, 'bold');
lbl(cg, 'Textbook rules applied before scoring.', 12, T.label2);
h.rules = lbl(cg, '', 12.5, T.label);
h.rules.WordWrap = 'on'; h.rules.VerticalAlignment = 'top';

% ================================================================ start
applyPreset('Balanced');
setSeg(h.propSeg, prop);
layoutOrbits();
highlightChip(1);
update(false);

% ================================================================ callbacks
    function loadExample(i)
        d = defaults;
        v = examples{i, 2};
        fn = fieldnames(v);
        for f = 1:numel(fn), d.(fn{f}) = v.(fn{f}); end
        h.alt1.Value = d.alt1; h.alt2.Value = d.alt2; h.di.Value = d.di_deg; h.dw.Value = d.dw_deg;
        orient = 'relative';
        if ~isempty(d.elements)
            orient = 'elements';
            for q = 1:6, h.el(q).Value = d.elements(q); end
        end
        h.apo1.Value = max(d.alt1, d.apo1); h.apo2.Value = max(d.alt2, d.apo2);     % max ignores NaN
        shape1 = 'circular'; if ~isnan(d.apo1), shape1 = 'elliptical'; end
        shape2 = 'circular'; if ~isnan(d.apo2), shape2 = 'elliptical'; end
        layoutOrbits();
        h.dvb.Value = d.dv_budget; h.mrg.Value = d.margin_pct; h.tof.Value = d.tof_max_h; h.rbm.Value = d.rb_max;
        h.ispc.Value = d.isp_chem; h.tw.Value = d.tw_chem; h.ispe.Value = d.isp_ep; h.acc.Value = d.accel_ep;
        prop = d.propulsion; setSeg(h.propSeg, prop);
        applyPreset(d.preset);
        manual = false;
        highlightChip(i);
        update(true);
    end

    function inputsEdited()
        highlightChip(0);
        update(false);
    end

    function pickProp(v)
        prop = v; setSeg(h.propSeg, v);
        highlightChip(0);
        update(false);
    end

    function pickOrient(v)
        orient = v;
        layoutOrbits();
        highlightChip(0);
        update(false);
    end

    function pickShape(k2, v)
        if k2 == 1, shape1 = v; else, shape2 = v; end
        layoutOrbits();
        highlightChip(0);
        update(false);
    end

    function layoutOrbits()
        % Rows of the Orbits card: orientation; shape, peri, apo, i, RAAN, argp (initial); the
        % same for the target; then dw and di (relative inputs) or the angle between the planes
        e1 = strcmp(shape1, 'elliptical'); e2 = strcmp(shape2, 'elliptical'); el = strcmp(orient, 'elements');
        rows = [1, 1 1 e1 el el el&&e1, 1 1 e2 el el el&&e2, ~el&&e2, ~el, el];
        h.orbGrid.RowHeight = num2cell(40*rows);
        rowObj = [h.apo1.Parent, h.el(1:3).Parent, h.apo2.Parent, h.el(4:6).Parent, h.dw.Parent, h.di.Parent, h.planes];
        set(rowObj, {'Visible'}, arrayfun(@onOff, logical(rows([4:7, 10:16])), 'UniformOutput', false).');
        left.RowHeight{2} = 8 + 40*sum(rows);
        setSeg(h.orientSeg, orient); setSeg(h.shape1, shape1); setSeg(h.shape2, shape2);
        if e1, h.alt1Lbl.Text = 'Initial periapsis'; else, h.alt1Lbl.Text = 'Initial altitude'; end
        if e2, h.alt2Lbl.Text = 'Target periapsis'; else, h.alt2Lbl.Text = 'Target altitude'; end
    end

    function pickPreset(v)
        applyPreset(v);
        highlightChip(0);
        update(false);
    end

    function applyPreset(name)
        i = find(strcmp(pnames, name), 1);
        for k2 = 1:4
            h.w(k2).Value = W(i, k2);
            h.wLbl(k2).Text = sprintf('%.2f', W(i, k2));
        end
        setSeg(h.presetSeg, name);
        h.customNote.Visible = 'off';
    end

    function weightMoving(k2, val)
        h.wLbl(k2).Text = sprintf('%.2f', val);
    end

    function weightSet()
        for k2 = 1:4, h.wLbl(k2).Text = sprintf('%.2f', h.w(k2).Value); end
        setSeg(h.presetSeg, '');
        h.customNote.Visible = 'on';
        highlightChip(0);
        update(false);
    end

    function rowPicked(e)
        if isempty(e.Indices) || isempty(order), return; end
        key = order{e.Indices(1, 1)};
        if strcmp(key, selected), return; end
        selected = key; manual = true;
        renderTable(); renderTrade(); renderGeo(true);
    end

    function highlightChip(i)
        for j = 1:numel(h.chips)
            if j == i
                set(h.chips(j), 'BackgroundColor', T.accent, 'FontColor', T.surface, 'FontWeight', 'bold');
            else
                set(h.chips(j), 'BackgroundColor', T.fill, 'FontColor', T.label, 'FontWeight', 'normal');
            end
        end
    end

% ================================================================ evaluation
    function update(animate)
        w = arrayfun(@(s) s.Value, h.w);
        if sum(w) <= 0, w = [1 1 1 1]; end
        apo1 = NaN; apo2 = NaN;          % an apoapsis below the periapsis is raised to it
        if strcmp(shape1, 'elliptical'), apo1 = max(h.apo1.Value, h.alt1.Value); h.apo1.Value = apo1; end
        if strcmp(shape2, 'elliptical'), apo2 = max(h.apo2.Value, h.alt2.Value); h.apo2.Value = apo2; end
        dw = 0; if strcmp(shape2, 'elliptical'), dw = h.dw.Value; end     % a circular target has no apse line
        el = [];
        if strcmp(orient, 'elements'), el = arrayfun(@(f) f.Value, h.el); end
        m = odst_mission('alt1', h.alt1.Value, 'alt2', h.alt2.Value, 'apo1', apo1, 'apo2', apo2, ...
            'dw_deg', dw, 'di_deg', h.di.Value, 'elements', el, ...
            'dv_budget', h.dvb.Value, 'margin_pct', h.mrg.Value, 'tof_max_h', h.tof.Value, 'propulsion', prop, ...
            'isp_chem', h.ispc.Value, 'tw_chem', h.tw.Value, 'isp_ep', h.ispe.Value, 'accel_ep', h.acc.Value, ...
            'rb_max', h.rbm.Value, 'weights', w);
        [opts, rec, adv] = odst_evaluate(m);
        h.planes.Text = sprintf('    Planes %.4g%s apart', m.di_deg, DEG);
        lastM = m; lastOpts = opts; lastRec = rec; lastAdv = adv;
        if ~manual || ~any(strcmp({opts.key}, selected))
            manual = false;
            if isempty(rec), selected = opts(1).key; else, selected = rec.key; end
        end
        renderVerdict(); renderTable(); renderTrade(); renderGeo(animate); renderRules();
    end

% ================================================================ renderers
    function renderVerdict()
        r = lastRec; m = lastM; opts = lastOpts;
        if isempty(r)
            h.vDot.FontColor = T.bad;
            h.vEyebrow.Text = 'No feasible manoeuvre';
            setChip([' ' WRN ' Every candidate breaks a limit '], T.bad, T.badSoft);
            h.vName.Text = 'Relax a constraint';
            h.vSub.Text = 'Raise the delta-v budget or the time limit, or add propulsion';
            for k2 = 1:5, h.stat(k2).Text = '-'; end
            msgs = {};
            for i = 1:numel(opts)
                if opts(i).applicable
                    msgs{end+1} = sprintf('%s: %s.', opts(i).name, opts(i).reasons); %#ok<AGROW>
                end
            end
            h.vWhy.Text = strjoin(msgs, '  ');
            return
        end
        f = find([opts.feasible]);
        [~, j] = max([opts(f).score_topsis]);
        topT = opts(f(j));
        h.vDot.FontColor = T.series.(r.key);
        h.vEyebrow.Text = 'Recommended manoeuvre';
        if strcmp(topT.key, r.key)
            setChip([' ' CHK ' TOPSIS agrees '], T.good, T.goodSoft);
            note = '';
        else
            setChip([' ' WRN ' Methods disagree '], T.warn, T.warnSoft);
            note = sprintf(' TOPSIS prefers %s (%.2f against %.2f), so this choice is close.', ...
                           topT.name, topT.score_topsis, r.score_topsis);
        end
        [a, b] = splitName(r.name);
        h.vName.Text = a; h.vSub.Text = b;
        vals = {sprintf('%.3f km/s', r.dv), odst_fmt_time(r.tof), burnsText(r.burns), ...
                sprintf('%.1f', 100*r.prop_frac), sprintf('%.2f', r.score_ws)};
        for k2 = 1:5, h.stat(k2).Text = vals{k2}; end
        h.vWhy.Text = [detailText(r, m) note sprintf(' %d of %d candidates meet your constraints.', ...
                       numel(f), numel(opts))];
    end

    function setChip(txt, fg, bg)
        set(h.vChip, 'Text', txt, 'FontColor', fg, 'BackgroundColor', bg);
    end

    function renderTable()
        opts = lastOpts;
        sc = [opts.score_ws]; sc(isnan(sc)) = -1;
        [~, idx] = sortrows([-double([opts.feasible]).', -sc.']);
        order = {opts(idx).key};
        D = cell(numel(idx), 9);
        for r = 1:numel(idx)
            o = opts(idx(r));
            if o.feasible
                score = sprintf('%.2f', o.score_ws); tops = sprintf('%.2f', o.score_topsis);
                if ~isempty(lastRec) && strcmp(o.key, lastRec.key)
                    status = [CHK ' Recommended'];
                else
                    status = 'Feasible';
                end
            else
                score = '-'; tops = '-'; status = o.reasons;
            end
            D(r, :) = {o.name, sprintf('%.3f', o.dv), odst_fmt_time(o.tof), burnsText(o.burns), ...
                       sprintf('%.1f%%', 100*o.prop_frac), sprintf('%.1f', o.complexity), score, tops, status};
            if isnan(o.dv), D(r, 2:5) = {'-'}; end      % an elliptical model that doesn't apply
        end
        h.tbl.Data = D;
        h.tbl.ColumnWidth = {250, 76, 76, 52, 80, 82, 58, 60, 'auto'};
        removeStyle(h.tbl);
        addStyle(h.tbl, uistyle('HorizontalAlignment', 'right'), 'column', 2:8);
        off = find(~[opts(idx).feasible]);
        if ~isempty(off)
            addStyle(h.tbl, uistyle('FontColor', T.label3), 'row', off);
            addStyle(h.tbl, uistyle('FontColor', T.bad), 'cell', [off(:), 9*ones(numel(off), 1)]);
        end
        sel = find(strcmp(order, selected));
        recRow = [];
        if ~isempty(lastRec), recRow = find(strcmp(order, lastRec.key)); end
        if ~isempty(sel) && ~isequal(sel, recRow)
            addStyle(h.tbl, uistyle('BackgroundColor', T.fill), 'row', sel);
        end
        if ~isempty(recRow)
            addStyle(h.tbl, uistyle('BackgroundColor', T.accentSoft, 'FontWeight', 'bold'), 'row', recRow);
        end
    end

    function renderTrade()
        ax = h.axTrade; opts = lastOpts; m = lastM;
        cla(ax); hold(ax, 'on');
        app = opts([opts.applicable]);
        x = max([app.tof]/3600, 1/60); y = [app.dv];
        xl = 10.^[floor(log10(min(x))) - 0.15, ceil(log10(max(x))) + 0.15];
        ytop = 1.15*max([y, m.dv_budget*(m.dv_budget < 1.6*max(y))]);
        for i = 1:numel(app)
            c = T.series.(app(i).key);
            sz = 70; if strcmp(app(i).key, selected), sz = 130; end
            if app(i).feasible
                scatter(ax, x(i), y(i), sz, c, 'filled', 'MarkerEdgeColor', T.surface, 'LineWidth', 1.5);
            else
                scatter(ax, x(i), y(i), sz, 'MarkerEdgeColor', c, 'MarkerFaceColor', T.surface, 'LineWidth', 1.8);
            end
        end
        % direct labels, nudged apart when two markers sit close together
        [~, o2] = sort(y, 'descend'); placed = zeros(0, 2);
        for i = o2
            lx = log10(x(i)); ly = y(i);
            for p = 1:size(placed, 1)
                if abs(placed(p, 1) - lx) < 0.9 && abs(placed(p, 2) - ly) < 0.06*ytop
                    ly = placed(p, 2) - 0.06*ytop;
                end
            end
            placed(end+1, :) = [lx, ly]; %#ok<AGROW>
            nearRight = lx > log10(xl(2)) - 0.9;
            if nearRight, tx = x(i)/1.3; ha = 'right'; else, tx = x(i)*1.3; ha = 'left'; end
            col = T.label; if ~app(i).feasible, col = T.label3; end
            text(ax, tx, ly, shortName(app(i).key), 'Color', col, 'FontSize', 11.5, ...
                 'FontName', T.font, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
        end
        if m.tof_max_h > xl(1) && m.tof_max_h < xl(2)
            xline(ax, m.tof_max_h, '--', 'Time limit', 'Color', T.bad, 'LineWidth', 1, 'FontSize', 10, ...
                  'LabelVerticalAlignment', 'bottom', 'FontName', T.font);
        end
        if m.dv_budget < ytop
            yline(ax, m.dv_budget, '--', [DEL 'v budget'], 'Color', T.bad, 'LineWidth', 1, 'FontSize', 10, ...
                  'FontName', T.font);
        end
        ax.XScale = 'log'; xlim(ax, xl); ylim(ax, [0 ytop]);
        ticks = [1/60 1/6 1 6 24 240 2400 8766 87660];
        names = {'1 min', '10 min', '1 h', '6 h', '1 day', '10 days', '100 days', '1 yr', '10 yr'};
        in = ticks >= xl(1) & ticks <= xl(2);
        ax.XTick = ticks(in); ax.XTickLabel = names(in);
        ax.XMinorGrid = 'off'; ax.XMinorTick = 'off';
        xlabel(ax, 'Transfer time (log scale)'); ylabel(ax, ['Total ' DEL 'v [km/s]']);
        hold(ax, 'off');
    end

    function renderGeo(animate)
        ax = h.axGeo; m = lastM; opts = lastOpts;
        o = opts(strcmp({opts.key}, selected));
        if isempty(o), o = opts(1); end
        h.export.Enable = onOff(o.applicable && ~isempty(o.vplan));     % low thrust has no burns to export
        r1 = C.RE + m.alt1; r2 = C.RE + m.alt2; k = 1e-3;
        ra1 = C.RE + max(m.alt1, m.apo1); ra2 = C.RE + max(m.alt2, m.apo2);   % max ignores NaN
        ell = any(strcmp({opts.key}, 'apse'));
        cla(ax); hold(ax, 'on');
        th = linspace(0, 2*pi, 361);
        fill(ax, k*C.RE*cos(th), k*C.RE*sin(th), T.earth, 'EdgeColor', 'none');
        if ell        % top view of both orbits; the target's plane is turned by di about x
            [x1, y1] = outline(r1, ra1, 0, m.w1_deg*pi/180); [x2, y2] = outline(r2, ra2, m.di_deg*pi/180, m.dw_deg*pi/180);
            plot(ax, k*x1, k*y1, 'Color', T.orbit, 'LineWidth', 1);
            plot(ax, k*x2, k*y2, '--', 'Color', T.orbit, 'LineWidth', 1);
            r1 = ra1; r2 = ra2;                                  % for the axis limits below
        else
            plot(ax, k*r1*cos(th), k*r1*sin(th), 'Color', T.orbit, 'LineWidth', 1);
            if abs(r2 - r1) > 1e-6
                plot(ax, k*r2*cos(th), k*r2*sin(th), '--', 'Color', T.orbit, 'LineWidth', 1);
            end
        end
        c = T.series.(o.key);
        [a, b] = splitName(o.name);
        if isempty(b), h.geoTitle.Text = a; else, h.geoTitle.Text = [a ', ' lower(b(1)) b(2:end)]; end
        if strcmp(o.key, 'low_thrust')
            f = linspace(0, 1, 700); rr = r1 + (r2 - r1)*f; aa = 12*pi*f;      % 6 turns drawn
            X = k*rr.*cos(aa); Y = k*rr.*sin(aa); bi = [];
            v1 = sqrt(C.MU/r1); v2 = sqrt(C.MU/r2);
            revs = abs(v1^4 - v2^4)/(8*pi*C.MU*m.accel_ep*1e-3);
            h.geoSub.Text = sprintf('Schematic spiral: 6 turns drawn, about %s in flight.', commas(revs));
        elseif ell && ~o.applicable
            X = []; Y = []; bi = [];
            h.geoSub.Text = ['Not available: ' o.note '.'];
        elseif ell
            [tt, YY] = odst_simulate_vec(o.y0, o.vplan, 0);
            X = k*YY(:, 1).'; Y = k*YY(:, 2).';
            bi = zeros(1, size(o.vplan, 1));
            for q = 1:size(o.vplan, 1), [~, bi(q)] = min(abs(tt - o.vplan(q, 1))); end
            h.geoSub.Text = 'Top view in 10^3 km. Solid: start orbit, dashed: target orbit.';
        elseif isempty(o.plan) || o.burns == 1
            X = k*r1; Y = 0; bi = 1;
            h.geoSub.Text = 'A single burn on the initial orbit rotates the plane.';
        else
            [tt, YY] = odst_simulate(r1, o.plan, 0);
            X = k*YY(:, 1).'; Y = k*YY(:, 2).';
            bi = zeros(1, size(o.plan, 1));
            for q = 1:size(o.plan, 1), [~, bi(q)] = min(abs(tt - o.plan(q, 1))); end
            h.geoSub.Text = sprintf('Top view in 10^3 km. Solid: %s km start orbit, dashed: %s km target.', ...
                                    commas(m.alt1), commas(m.alt2));
        end
        lim = 1.12*max([k*r2, k*r1, abs(X), abs(Y)]);
        axis(ax, 'equal'); xlim(ax, [-lim lim]); ylim(ax, [-lim lim]);
        xlabel(ax, 'x [10^3 km]'); ylabel(ax, 'y [10^3 km]');
        ln = plot(ax, NaN, NaN, 'Color', c, 'LineWidth', 2.2);
        mk = gobjects(0); tx = gobjects(0);
        for q = 1:numel(bi)
            mk(q) = plot(ax, X(bi(q)), Y(bi(q)), 'o', 'MarkerSize', 8, 'MarkerFaceColor', c, ...
                         'MarkerEdgeColor', T.surface, 'LineWidth', 1.5, 'Visible', 'off');
            tx(q) = text(ax, X(bi(q)), Y(bi(q)), sprintf('  %sv%d', DEL, q), 'Color', T.label2, ...
                         'FontSize', 11, 'FontWeight', 'bold', 'FontName', T.font, ...
                         'VerticalAlignment', 'bottom', 'Visible', 'off');
        end
        if numel(X) > 1
            s = [0, cumsum(hypot(diff(X), diff(Y)))];
            frac = s(bi)/s(end);
            if animate
                % progress follows a critically damped spring (response 1.05 s),
                % the same timing as the web planner's draw-in
                craft = plot(ax, X(1), Y(1), 'o', 'MarkerSize', 7, 'MarkerFaceColor', T.surface, ...
                             'MarkerEdgeColor', c, 'LineWidth', 2);
                wn = 2*pi/1.05; t0 = tic; p = 0;
                while p < 0.995 && isvalid(ln)
                    t = toc(t0);
                    p = 1 - (1 + wn*t)*exp(-wn*t);
                    j = find(s <= p*s(end), 1, 'last');
                    set(ln, 'XData', X(1:j), 'YData', Y(1:j));
                    set(craft, 'XData', X(j), 'YData', Y(j));
                    on = frac <= p + 1e-3;
                    set(mk(on), 'Visible', 'on'); set(tx(on), 'Visible', 'on');
                    drawnow limitrate
                end
                if isvalid(craft), delete(craft); end
            end
            if isvalid(ln), set(ln, 'XData', X, 'YData', Y); end
        end
        set(mk(isvalid(mk)), 'Visible', 'on'); set(tx(isvalid(tx)), 'Visible', 'on');
        if isvalid(ax), hold(ax, 'off'); end
    end

    function exportPlan()
        % the selected transfer as a GMAT script plus a JSON file of its burns, same base name
        o = lastOpts(strcmp({lastOpts.key}, selected));
        [f, p] = uiputfile('*.script', 'Export the transfer to NASA GMAT', ['odst_' o.key '.script']);
        if isequal(f, 0), return; end
        [~, base] = fileparts(f);
        writeText(fullfile(p, f), odst_gmat_script(lastM, o));
        writeText(fullfile(p, [base '.json']), odst_burns_json(lastM, o));
        h.geoSub.Text = sprintf('Saved %s and %s.json.', f, base);
    end

    function writeText(file, txt)
        fid = fopen(file, 'w', 'n', 'UTF-8');
        fprintf(fid, '%s', txt);
        fclose(fid);
    end

    function [x, y] = outline(rp, ra, di, dw)
        ang = linspace(0, 2*pi, 361);
        x = zeros(size(ang)); y = x;
        for q = 1:numel(ang)
            r = odst_orbit_state(rp, ra, di, dw, ang(q));
            x(q) = r(1); y(q) = r(2);
        end
    end

    function renderRules()
        h.rules.Text = strcat({[BUL '  ']}, lastAdv);
    end

% ================================================================ text helpers
    function s = detailText(o, m)
        p = o.params;
        switch o.key
            case 'hohmann_split'
                s = sprintf('Split the %g%s plane change: %.2f%s at departure and %.2f%s at arrival.', ...
                    m.di_deg, DEG, p.split_at_r1*m.di_deg, DEG, (1 - p.split_at_r1)*m.di_deg, DEG);
            case 'bielliptic'
                s = sprintf(['Coast out to an optimal apoapsis of %s km (%.1f times the initial radius), ' ...
                    'where a plane change is cheapest, then fall to the target orbit.'], commas(p.rb), p.rb_over_r1);
            case 'one_tangent'
                s = sprintf('A faster transfer ellipse (k = %.4f, e = %.3f) meets the target orbit at a %.1f%s flight-path angle.', ...
                    p.k, p.e, p.fpa_deg, DEG);
                if startsWith(o.note, 'k solved')
                    s = [s ' k is the smallest value that meets the time limit.'];
                end
            case 'low_thrust'
                s = sprintf('Continuous thrust at %g m/s%s spirals out to the target. The high Isp keeps the propellant mass low.', ...
                    m.accel_ep, char(178));
            case 'apse'
                if o.burns == 1
                    s = 'One burn where the two orbits touch puts the spacecraft on the target orbit.';
                    if m.di_deg > 0
                        s = sprintf('%s It also turns the plane by %g%s.', s, m.di_deg, DEG);
                    end
                else
                    s = sprintf('Two burns half a transfer ellipse apart, at %s km and %s km from Earth''s centre.', ...
                        commas(p.rA), commas(p.rB));
                    if m.di_deg > 0
                        s = sprintf('%s %.2f%s of the %g%s plane change is done at departure.', s, ...
                            p.split_at_dep*m.di_deg, DEG, m.di_deg, DEG);
                    end
                end
            case 'lambert'
                s = sprintf(['The cheapest two-burn transfer the Lambert search found: depart at %.0f%s true ' ...
                    'anomaly on the start orbit, arrive at %.0f%s on the target orbit.'], p.th1_deg, DEG, p.th2_deg, DEG);
                if p.time_limited, s = [s ' The time limit sets the transfer time.']; end
            otherwise
                if o.burns == 1
                    s = 'A single burn at the node rotates the orbit plane. No transfer orbit is needed.';
                elseif m.di_deg > 0
                    s = sprintf('The %g%s plane change is combined with the burn at the larger radius, where the orbital speed is lowest.', ...
                        m.di_deg, DEG);
                else
                    s = 'Two tangential burns at the ends of the transfer ellipse.';
                end
        end
    end

    function n = shortName(key)
        switch key
            case 'hohmann',       n = 'Hohmann';
            case 'hohmann_split', n = 'Hohmann split';
            case 'bielliptic',    n = 'Bi-elliptic';
            case 'one_tangent',   n = 'One-tangent';
            case 'apse',          n = 'Apse-to-apse';
            case 'lambert',       n = 'Lambert';
            otherwise,            n = 'Low-thrust';
        end
    end

    function t = burnsText(b)
        if b == 0, t = 'Cont.'; else, t = sprintf('%d', b); end
    end

    function [a, b] = splitName(n)
        q = strfind(n, ' (');
        if isempty(q)
            a = n; b = '';
        else
            a = n(1:q(1) - 1); b = n(q(1) + 2:end - 1);
            b(1) = upper(b(1));
            b = strrep(strrep(b, ' di ', [' ' DEL 'i ']), ' dv', [' ' DEL 'v']);
        end
    end

    function s = commas(x)
        s = sprintf('%d', round(x));
        for q = numel(s) - 3:-3:1
            s = [s(1:q) ',' s(q + 1:end)];
        end
    end

% ================================================================ UI builders
    function L = lbl(parent, txt, fs, color, weight)
        if nargin < 5, weight = 'normal'; end
        L = uilabel(parent, 'Text', txt, 'FontName', T.font, 'FontSize', fs, ...
                    'FontColor', color, 'FontWeight', weight);
    end

    function groupTitle(parent, txt)
        lbl(parent, ['  ' txt], 12.5, T.label2, 'bold');
    end

    function g2 = listPanel(parent, nrows)
        p = uipanel(parent, 'BorderType', 'none', 'BackgroundColor', T.surface);
        g2 = uigridlayout(p, [nrows 1], 'RowHeight', repmat({40}, 1, nrows), 'Padding', [14 4 10 4], ...
                          'RowSpacing', 0, 'BackgroundColor', T.surface);
    end

    function [f, L] = numRow(g2, label, value, limits, unit)
        r = uigridlayout(g2, [1 3], 'ColumnWidth', {'1x', 96, 42}, 'Padding', [0 0 0 0], ...
                         'ColumnSpacing', 6, 'BackgroundColor', T.surface);
        L = lbl(r, label, 13.5, T.label);
        f = uieditfield(r, 'numeric', 'Value', value, 'Limits', limits, 'HorizontalAlignment', 'right', ...
                        'ValueDisplayFormat', '%.10g', ...
                        'FontName', T.font, 'FontSize', 13.5, 'FontColor', T.label, ...
                        'BackgroundColor', T.surface, 'ValueChangedFcn', @(~, ~) inputsEdited());
        lbl(r, unit, 12, T.label3);
    end

    function seg = segRow(g2, label, labels, values, onPick)
        r = uigridlayout(g2, [1 2], 'ColumnWidth', {'1x', 170}, 'Padding', [0 4 0 4], ...
                         'ColumnSpacing', 6, 'BackgroundColor', T.surface);
        lbl(r, label, 13.5, T.label);
        seg = segmented(r, labels, values, onPick);
    end

    function seg = shapeRow(g2, label, k2)
        r = uigridlayout(g2, [1 2], 'ColumnWidth', {'1x', 170}, 'Padding', [0 4 0 4], ...
                         'ColumnSpacing', 6, 'BackgroundColor', T.surface);
        lbl(r, label, 13.5, T.label);
        seg = segmented(r, {'Circular', 'Elliptical'}, {'circular', 'elliptical'}, @(v) pickShape(k2, v));
    end

    function s = onOff(tf)
        if tf, s = 'on'; else, s = 'off'; end
    end

    function [s, valueLbl] = sliderRow(g2, label, value, k2)
        r = uigridlayout(g2, [1 3], 'ColumnWidth', {92, '1x', 36}, 'Padding', [4 0 4 0], ...
                         'ColumnSpacing', 8, 'BackgroundColor', T.surface);
        lbl(r, label, 13, T.label);
        s = uislider(r, 'Limits', [0 1], 'Value', value, 'MajorTicks', [], 'MinorTicks', [], ...
                     'ValueChangingFcn', @(~, e) weightMoving(k2, e.Value), ...
                     'ValueChangedFcn', @(~, ~) weightSet());
        valueLbl = lbl(r, sprintf('%.2f', value), 12.5, T.label2);
        valueLbl.HorizontalAlignment = 'right';
    end

    function seg = segmented(parent, labels, values, onPick)
        % iOS-style segmented control: grey track, selected segment white and bold
        tr = uigridlayout(parent, [1 numel(labels)], 'Padding', [2 2 2 2], 'ColumnSpacing', 2, ...
                          'RowHeight', {'1x'}, 'BackgroundColor', T.fill);
        seg = struct('values', {values}, 'buttons', gobjects(1, numel(labels)));
        for q = 1:numel(labels)
            v = values{q};
            seg.buttons(q) = uibutton(tr, 'Text', labels{q}, 'FontName', T.font, 'FontSize', 12.5, ...
                'FontColor', T.label, 'BackgroundColor', T.fill, 'ButtonPushedFcn', @(~, ~) onPick(v));
        end
    end

    function setSeg(seg, value)
        for q = 1:numel(seg.buttons)
            if strcmp(seg.values{q}, value)
                set(seg.buttons(q), 'BackgroundColor', T.surface, 'FontWeight', 'bold');
            else
                set(seg.buttons(q), 'BackgroundColor', T.fill, 'FontWeight', 'normal');
            end
        end
    end

    function [ax, sub, ttl, btn] = chartCard(parent, ttlText, subText)
        p = uipanel(parent, 'BorderType', 'none', 'BackgroundColor', T.surface);
        g2 = uigridlayout(p, [3 1], 'RowHeight', {22, 18, '1x'}, 'Padding', [18 12 18 10], ...
                          'RowSpacing', 2, 'BackgroundColor', T.surface);
        hd = uigridlayout(g2, [1 2], 'ColumnWidth', {'1x', 'fit'}, 'Padding', [0 0 0 0], ...
                          'ColumnSpacing', 8, 'BackgroundColor', T.surface);
        ttl = lbl(hd, ttlText, 16, T.label, 'bold');
        if nargout > 3          % a small button at the right of the title
            btn = uibutton(hd, 'FontName', T.font, 'FontSize', 12, 'BackgroundColor', T.fill, ...
                           'FontColor', T.label);
        end
        sub = lbl(g2, subText, 12, T.label2);
        ax = uiaxes(g2);
        set(ax, 'Color', T.surface, 'BackgroundColor', T.surface, 'FontName', T.font, 'FontSize', 11, ...
            'XColor', T.label2, 'YColor', T.label2, 'Box', 'off', 'TickLength', [0 0], 'LineWidth', 0.8);
        set(ax, 'GridColor', T.grid, 'GridAlpha', 1, 'XGrid', 'on', 'YGrid', 'on');
        set([ax.XLabel, ax.YLabel], 'Color', T.label2);
        try
            ax.Toolbar.Visible = 'off';
            disableDefaultInteractivity(ax);
        catch
        end
    end
end
