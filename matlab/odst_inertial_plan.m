function [y0, vplan] = odst_inertial_plan(m, o)
%ODST_INERTIAL_PLAN Option o's start state and vector burn plan in EarthMJ2000Eq.
%   [y0, vplan] = odst_inertial_plan(m, o): the relative frame turned by m.R for a
%   mission given by orbital elements (with J2, the frame of the mission o was planned
%   for, o.m_eff); unchanged otherwise (the tool's frame). Mirrors inertial_plan() in
%   python/odst.py.
y0 = o.y0; vplan = o.vplan; R = m.R;
if isfield(o, 'm_eff') && ~isempty(o.m_eff), R = o.m_eff.R; end
if isempty(R), return; end
y0 = [R*y0(1:3); R*y0(4:6)];
vplan = [vplan(:, 1), (R*vplan(:, 2:4).').'];
end
