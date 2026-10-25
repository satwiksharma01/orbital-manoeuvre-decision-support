function [y0, vplan] = odst_inertial_plan(m, o)
%ODST_INERTIAL_PLAN Option o's start state and vector burn plan in EarthMJ2000Eq.
%   [y0, vplan] = odst_inertial_plan(m, o): the relative frame turned by m.R for a
%   mission given by orbital elements; unchanged otherwise (the tool's frame).
%   Mirrors inertial_plan() in python/odst.py.
y0 = o.y0; vplan = o.vplan;
if isempty(m.R), return; end
y0 = [m.R*y0(1:3); m.R*y0(4:6)];
vplan = [vplan(:, 1), (m.R*vplan(:, 2:4).').'];
end
