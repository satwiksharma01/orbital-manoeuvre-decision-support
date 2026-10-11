function style_figure(fig)
%STYLE_FIGURE Apply the ODST visual style to a figure before export.
%   style_figure(fig) uses the same design tokens as the web planner and the
%   GUI (odst_theme): white canvas, near-black titles aligned left, grey axis
%   text, hairline grid, no box or tick marks, no minor grid. It forces the
%   light theme where MATLAB supports themes (R2025a+ follows a dark desktop
%   otherwise). Series colours are left unchanged. Call it just before
%   exportgraphics.
if nargin < 1 || isempty(fig)
    fig = gcf;
end
T = odst_theme();
if isprop(fig, 'Theme')
    fig.Theme = 'light';
end
fig.Color = T.surface;

% Titles of tiledlayout figures (e.g. the Case A/B/C bar charts)
for tl = findall(fig, 'Type', 'tiledlayout').'
    set(tl.Title, 'Color', T.label, 'FontName', T.font, 'FontSize', 16, 'FontWeight', 'bold');
    if isprop(tl, 'Subtitle')
        set(tl.Subtitle, 'Color', T.label2, 'FontName', T.font, 'FontSize', 12, 'FontWeight', 'normal');
    end
end

for ax = findall(fig, 'Type', 'axes').'
    % axis colour first: grid colours are set afterwards so they stay hairline grey
    set(ax, 'Color', T.surface, 'FontName', T.font, 'FontSize', 11, ...
        'XColor', T.label2, 'YColor', T.label2, 'ZColor', T.label2, ...
        'Box', 'off', 'TickLength', [0 0], 'LineWidth', 0.8, 'Layer', 'bottom');
    set(ax, 'GridColor', T.grid, 'GridAlpha', 1, ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off', 'XMinorTick', 'off', 'YMinorTick', 'off');
    if ~isempty(findobj(ax, 'Type', 'bar'))
        ax.XGrid = 'off';                       % categories need no vertical lines
    end
    try
        ax.YAxis(1).Axle.Visible = 'off';       % baseline only, like the web charts
    catch
    end

    % text annotations (findobj skips the hidden title / label handles)
    set(findobj(ax, 'Type', 'text'), 'Color', T.label, 'FontName', T.font);

    set(ax.Title, 'Color', T.label, 'FontName', T.font, 'FontSize', 13, 'FontWeight', 'bold');
    if isprop(ax, 'TitleHorizontalAlignment')
        ax.TitleHorizontalAlignment = 'left';
    end
    if isprop(ax, 'Subtitle')
        set(ax.Subtitle, 'Color', T.label2, 'FontName', T.font, 'FontWeight', 'normal');
    end
    set([ax.XLabel, ax.YLabel], 'Color', T.label2, 'FontName', T.font, 'FontSize', 11);

    lg = ax.Legend;
    if ~isempty(lg)
        set(lg, 'TextColor', T.label, 'Color', T.surface, 'Box', 'off', ...
            'FontName', T.font, 'FontSize', 10.5);
    end
end

% Reference lines (xline / yline) and their labels
set(findall(fig, 'Type', 'constantline'), 'Color', T.label2);

for cb = findall(fig, 'Type', 'colorbar').'
    set(cb, 'Color', T.label2, 'TickLength', 0, 'Box', 'off', 'FontName', T.font, 'FontSize', 11);
    cb.Label.Color = T.label2;
end
end
