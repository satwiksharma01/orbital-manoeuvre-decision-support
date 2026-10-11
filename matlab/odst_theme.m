function T = odst_theme()
%ODST_THEME Shared ODST design tokens (same palette as web/index.html).
%   T = odst_theme() returns the colours (RGB triplets in [0 1]) and the font
%   used by the GUI (odst_app) and every exported figure (style_figure), so
%   the MATLAB outputs match the web planner and the report.
%
%   T.bg, T.surface          page and card backgrounds
%   T.label, T.label2/3      primary / secondary / tertiary text
%   T.grid, T.sep, T.fill    hairlines and control fills
%   T.accent, T.good, ...    semantic colours (+ *Soft tints for chips)
%   T.series.<key>           fixed colour per manoeuvre (never by rank)
%   T.font                   first installed face of an Apple-like stack
persistent cached
if ~isempty(cached)
    T = cached;
    return
end
hex = @(h) sscanf(h(2:end), '%2x%2x%2x', [1 3])/255;

T.bg         = hex('#f5f5f7');
T.surface    = hex('#ffffff');
T.label      = hex('#1d1d1f');
T.label2     = hex('#6e6e73');
T.label3     = hex('#86868b');
T.sep        = hex('#e3e3e6');
T.grid       = hex('#ececf0');
T.fill       = hex('#ececee');
T.accent     = hex('#0071e3');
T.accentSoft = hex('#e8f2fd');
T.good       = hex('#248a3d');
T.goodSoft   = hex('#e7f2ea');
T.warn       = hex('#b25000');
T.warnSoft   = hex('#fff0db');
T.bad        = hex('#d70015');
T.badSoft    = hex('#fcebec');
T.earth      = hex('#a9cbf3');
T.orbit      = hex('#aeaeb2');
T.none       = hex('#e5e5ea');
T.series = struct('hohmann', hex('#2a78d6'), 'hohmann_split', hex('#eb6834'), ...
                  'bielliptic', hex('#1baf7a'), 'one_tangent', hex('#eda100'), ...
                  'low_thrust', hex('#e87ba4'));
% Ordinal ramp of the bi-elliptic hue (light to dark) for its rB variants
T.bielliptic_ramp = [hex('#46c094'); hex('#1baf7a'); hex('#13875e'); hex('#0b5c40')];
T.font = pickFont({'SF Pro Text', 'Helvetica Neue', 'Segoe UI', 'Inter', 'Helvetica', ...
                   'Arial', 'Liberation Sans', 'Noto Sans', 'DejaVu Sans'});
cached = T;
end

function f = pickFont(stack)
% First font of the stack that is installed; MATLAB's default otherwise.
f = get(groot, 'defaultAxesFontName');
try
    have = listfonts;
catch
    return
end
for i = 1:numel(stack)
    if any(strcmpi(have, stack{i}))
        f = stack{i};
        return
    end
end
end
