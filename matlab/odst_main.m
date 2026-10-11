%% ODST_MAIN  Orbital Manoeuvre Decision Support Tool - full analysis run
%  ASE3004 Space Flight Mechanics, Project Option 38.
%
%  Runs: (1) validation suite, (2) three mission case studies,
%        (3) analysis figures, (4) Simulink finite-burn verification
%        (if Simulink is installed). Interactive version: odst_app
%
%  Figures use the shared ODST design (odst_theme, style_figure) and are
%  exported at 300 DPI to matlab/results.
%
%  Requires MATLAB R2020b or newer (base MATLAB only; Simulink optional).
clear; clc; close all;
C = odst_const(); RE = C.RE;
TH = odst_theme();
[pnames, W] = odst_presets();
outDir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% 1. Validation
nFail = odst_tests();

%% 2. Mission case studies
caseA = odst_mission('alt1', 400, 'alt2', 35786, 'di_deg', 28.5, 'dv_budget', 6, ...
                     'tof_max_h', 24*200, 'weights', W(1, :));                 % Balanced
caseB = odst_mission('alt1', 7000 - RE, 'alt2', 105000 - RE, 'di_deg', 0, 'dv_budget', 6, ...
                     'tof_max_h', 24*60, 'propulsion', 'chemical', 'weights', W(2, :));
caseC = odst_mission('alt1', 400, 'alt2', 800, 'di_deg', 0, 'dv_budget', 1, ...
                     'tof_max_h', 0.6, 'weights', W(3, :));                    % Time-critical
cases = {caseA, caseB, caseC};
titles = {'Case A: LEO 400 km \rightarrow GEO, \Deltai = 28.5\circ (Balanced)', ...
          'Case B: 7000 \rightarrow 105 000 km, R = 15 (Propellant-critical, chemical only)', ...
          'Case C: LEO 400 \rightarrow 800 km, TOF \leq 36 min (Time-critical)'};
caseNames = {'Case A', 'Case B', 'Case C'};
results = cell(size(cases));
for c = 1:numel(cases)
    [opts, rec, adv] = odst_evaluate(cases{c});
    fprintf('\n##### %s', caseNames{c});
    odst_print(cases{c}, opts, rec, adv);
    results{c} = opts;
    plotCase(opts, rec, titles{c});
    saveFig(gcf, fullfile(outDir, sprintf('case_%c_bars.png', 'A' + c - 1)));
end

%% 3a. Propagated trajectories of every impulsive option, case B
reportFigure('Name', 'Trajectories', 'Position', [60 60 1300 430]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact');
title(tl, 'Numerically propagated transfers (two-body, impulsive): Case B');
show = {'hohmann', 'bielliptic', 'one_tangent'};
for s = 1:numel(show)
    o = results{2}(strcmp({results{2}.key}, show{s}));
    nexttile; plotTransfer(RE + caseB.alt1, RE + caseB.alt2, o);
end
saveFig(gcf, fullfile(outDir, 'case_B_trajectories.png'));

%% 3b. Rule R1: Hohmann vs bi-elliptic as a function of R
r1 = 7000; R = logspace(log10(1.5), log10(80), 300);
reportFigure('Name', 'Hohmann vs bi-elliptic'); hold on; grid on;
plot(R, arrayfun(@(x) odst_hohmann(r1, x*r1, 0), R)/sqrt(C.MU/r1), 'LineWidth', 2, ...
     'Color', palette('hohmann'), 'DisplayName', 'Hohmann');
rbs = [2 5 15 1e5];
for n = 1:numel(rbs)
    db = arrayfun(@(x) odst_bielliptic(r1, x*r1, rbs(n)*x*r1, 0), R)/sqrt(C.MU/r1);
    if rbs(n) > 1e4, lab = 'Bi-elliptic, r_B = \infty'; ls = '--'; else, lab = sprintf('Bi-elliptic, r_B = %g r_2', rbs(n)); ls = '-'; end
    plot(R, db, ls, 'LineWidth', 1.8, 'Color', TH.bielliptic_ramp(n, :), 'DisplayName', lab);
end
xline(11.94, ':', 'R = 11.94', 'HandleVisibility', 'off');
xline(15.58, ':', 'R = 15.58', 'HandleVisibility', 'off');
set(gca, 'XScale', 'log'); ylim([0.4 0.58]);
xlabel('Radius ratio R = r_2/r_1'); ylabel('\Deltav_{total} / v_{c1}');
title('Hohmann vs bi-elliptic'); subtitle('Coplanar transfers, total \Deltav normalised by the initial circular speed');
legend('Location', 'southeast');
saveFig(gcf, fullfile(outDir, 'hohmann_vs_bielliptic.png'));

%% 3c. Numerical optimisation of the plane-change split (case A)
r1 = RE + caseA.alt1; r2 = RE + caseA.alt2; di = caseA.di_deg*pi/180;
s = linspace(0, 1, 401);
dvs = arrayfun(@(x) odst_hohmann(r1, r2, di, x), s);
sOpt = odst_hohmann_split_opt(r1, r2, di);
dvOpt = odst_hohmann(r1, r2, di, sOpt);
reportFigure('Name', 'Plane change split'); hold on; grid on;
plot(s*caseA.di_deg, dvs, 'LineWidth', 2, 'Color', palette('hohmann_split'));
plot(sOpt*caseA.di_deg, dvOpt, 'o', 'MarkerSize', 9, 'MarkerFaceColor', palette('hohmann_split'), ...
     'Color', TH.surface, 'LineWidth', 1.5);
text(1, max(dvs), {sprintf('optimum (marker): %.2f\\circ at the LEO burn', sOpt*caseA.di_deg), ...
     sprintf('\\Deltav = %.4f km/s, saves %.0f m/s', dvOpt, 1000*(dvs(1) - dvOpt))}, ...
     'VerticalAlignment', 'top');
xlabel('Plane change performed at the LEO (perigee) burn [\circ]'); ylabel('Total \Deltav [km/s]');
title(sprintf('Optimal plane-change split (\\Deltai = %.1f\\circ)', caseA.di_deg));
subtitle('Case A, Hohmann transfer LEO \rightarrow GEO');
saveFig(gcf, fullfile(outDir, 'plane_change_split.png'));

%% 3d. Decision map: recommended manoeuvre over (R, TOF limit)
keys = {'hohmann', 'hohmann_split', 'bielliptic', 'one_tangent', 'low_thrust'};
Rs = logspace(log10(1.05), log10(50), 30);
tofs = logspace(log10(0.3), log10(24*365), 26);
Z = nan(numel(tofs), numel(Rs));
for i = 1:numel(tofs)
    for j = 1:numel(Rs)
        m = odst_mission('alt1', 400, 'alt2', Rs(j)*(RE + 400) - RE, 'di_deg', 28.5, ...
                         'dv_budget', 10, 'tof_max_h', tofs(i), 'propulsion', 'both', 'weights', W(2, :));
        [~, rec] = odst_evaluate(m);
        if isempty(rec), Z(i, j) = 6; else, Z(i, j) = find(strcmp(keys, rec.key)); end
    end
end
reportFigure('Name', 'Decision map', 'Position', [100 100 900 620]);
imagesc(log10(Rs), log10(tofs), Z); axis xy;
cmap = [cell2mat(cellfun(@palette, keys, 'UniformOutput', false).'); TH.none];
colormap(cmap); caxis([0.5 6.5]); %#ok<CAXIS>
cb = colorbar('Ticks', 1:6, 'TickLabels', {'Hohmann', 'H. split', 'Bi-elliptic', 'One-tangent', 'Low-thrust', 'None'});
cb.Label.String = 'Recommended manoeuvre';
xlabel('log_{10} R = log_{10}(r_2/r_1)'); ylabel('log_{10} TOF limit [h]');
title('Decision map: recommended manoeuvre');
subtitle('\Deltai = 28.5\circ, propellant-critical weights, r_1 = LEO 400 km, chemical + electric');
saveFig(gcf, fullfile(outDir, 'decision_map.png'));

%% 3e. Weight sensitivity (case A)
wt = linspace(0, 0.9, 91);
base = caseA.weights; others = base([1 2 4])/sum(base([1 2 4]));
S = nan(numel(wt), numel(results{1}));
for i = 1:numel(wt)
    m = caseA; m.weights = [others(1:2)*(1 - wt(i)), wt(i), others(3)*(1 - wt(i))];
    o = odst_evaluate(m);
    for k = 1:numel(o)
        if o(k).feasible, S(i, k) = o(k).score_ws; end
    end
end
reportFigure('Name', 'Weight sensitivity'); hold on; grid on;
for k = 1:numel(results{1})
    if any(~isnan(S(:, k)))
        plot(wt, S(:, k), 'LineWidth', 2, 'Color', palette(results{1}(k).key), 'DisplayName', results{1}(k).name);
    end
end
xlabel('Weight on transfer time (others scaled proportionally)'); ylabel('Weighted-sum score');
title('Sensitivity to the transfer-time weight');
subtitle('Case A, LEO \rightarrow GEO: one-tangent overtakes Hohmann above a weight of about 0.86');
legend('Location', 'southoutside');
saveFig(gcf, fullfile(outDir, 'weight_sensitivity.png'));

%% 4. Simulink verification with finite burns (optional)
if license('test', 'Simulink') && ~isempty(ver('simulink'))
    [~, recA] = odst_evaluate(caseA);
    odst_run_simulink(recA, RE + caseA.alt1, RE + caseA.alt2, caseA.di_deg, 60);
    saveFig(gcf, fullfile(outDir, 'simulink_verification.png'));
else
    fprintf('\nSimulink not available - skipping finite-burn verification.\n');
end
fprintf('\nDone. Figures written to %s\n', outDir);

%% ------------------------------------------------------------ local functions
function f = reportFigure(varargin)
% White figure in the light theme (MATLAB R2025a+ otherwise follows a dark desktop).
f = figure(varargin{:});
if isprop(f, 'Theme'), f.Theme = 'light'; end
f.Color = 'w';
end

function saveFig(fig, file)
% Shared ODST style (odst_theme) + 300 DPI export.
style_figure(fig);
drawnow;
exportgraphics(fig, file, 'Resolution', 300, 'BackgroundColor', 'white');
end

function c = palette(key)
% Fixed colour per manoeuvre (follows the entity, never its rank).
T = odst_theme();
if isfield(T.series, key), c = T.series.(key); else, c = T.label3; end
end

function plotCase(opts, rec, ttl)
reportFigure('Name', ttl, 'Position', [60 60 1400 360]);
tl = tiledlayout(1, 4, 'TileSpacing', 'compact');
title(tl, ttl);
labels = {'Hohmann', 'H. split', 'Bi-elliptic', 'One-tangent', 'Low-thrust'};
keys = {opts.key};
lab = cellfun(@(k) labels{strcmp({'hohmann', 'hohmann_split', 'bielliptic', 'one_tangent', 'low_thrust'}, k)}, ...
              keys, 'UniformOutput', false);
sc = [opts.score_ws]; sc(~[opts.feasible]) = 0;
data = {[opts.dv], [opts.prop_frac], [opts.tof]/3600, sc};
names = {'Total \Deltav [km/s]', 'Propellant fraction m_p/m_0', 'Transfer time [h]', 'Weighted-sum score'};
for p = 1:4
    nexttile; hold on; grid on;
    b = bar(data{p}, 0.62, 'FaceColor', 'flat', 'EdgeColor', 'w', 'LineWidth', 1.5);
    for k = 1:numel(opts)
        b.CData(k, :) = palette(keys{k});
        if ~opts(k).feasible, b.CData(k, :) = 0.75 + 0.25*b.CData(k, :); end
    end
    if p == 3, set(gca, 'YScale', 'log'); end
    set(gca, 'XTick', 1:numel(opts), 'XTickLabel', lab, 'XTickLabelRotation', 30);
    title(names{p});
    if p == 4 && ~isempty(rec)
        k = find(strcmp(keys, rec.key));
        text(k, rec.score_ws + 0.05, 'recommended', 'HorizontalAlignment', 'center');
        ylim([0 1.15]);
    end
end
end

function plotTransfer(r1, r2, o)
C = odst_const(); T = odst_theme();
[tt, Y] = odst_simulate(r1, o.plan, 0);
th = linspace(0, 2*pi, 400);
hold on; axis equal; grid on;
k = 1e-3;                                   % plot in 10^3 km (no axis exponent)
fill(k*C.RE*cos(th), k*C.RE*sin(th), T.earth, 'EdgeColor', 'none');
plot(k*r1*cos(th), k*r1*sin(th), 'Color', T.orbit, 'LineWidth', 1);
plot(k*r2*cos(th), k*r2*sin(th), '--', 'Color', T.orbit, 'LineWidth', 1);
plot(k*Y(:, 1), k*Y(:, 2), 'LineWidth', 2.2, 'Color', palette(o.key));
for b = 1:size(o.plan, 1)                   % burn markers
    [~, i] = min(abs(tt - o.plan(b, 1)));
    plot(k*Y(i, 1), k*Y(i, 2), 'o', 'MarkerSize', 7, 'MarkerFaceColor', palette(o.key), ...
         'MarkerEdgeColor', T.surface, 'LineWidth', 1.2);
end
title(o.name);
subtitle(sprintf('\\Deltav = %.3f km/s, TOF = %s', o.dv, odst_fmt_time(o.tof)));
xlabel('x [10^3 km]'); ylabel('y [10^3 km]');
end
