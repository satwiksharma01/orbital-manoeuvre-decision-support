function m = odst_mission(varargin)
%ODST_MISSION Default mission definition, optionally overridden by name/value pairs.
%   m = odst_mission('alt1', 400, 'alt2', 35786, 'di_deg', 28.5)
%
%   alt1, alt2     initial / target altitude, the periapsis if elliptical [km]
%   apo1, apo2     initial / target apoapsis altitude [km]; NaN = circular
%   dw_deg         target periapsis, from the node line (the start periapsis if w1_deg = 0) [deg]
%   di_deg         inclination change: the angle between the orbit planes [deg]
%   w1_deg         start periapsis, from the node line [deg]
%   elements       [i1 RAAN1 argp1 i2 RAAN2 argp2] in deg, or [] (default). If set,
%                  di_deg, w1_deg and dw_deg follow from them, and R holds the rotation
%                  from the relative frame to EarthMJ2000Eq (odst_relative_geometry)
%   dv_budget      available delta-v [km/s]
%   margin_pct     delta-v margin the budget must also cover [%]
%   tof_max_h      maximum transfer time [h]
%   propulsion     'both' | 'chemical' | 'electric'
%   isp_chem       chemical specific impulse [s]
%   tw_chem        chemical thrust-to-weight at departure, T/(m0*g0)
%   isp_ep         electric specific impulse [s]
%   accel_ep       low-thrust acceleration F/m [m/s^2]
%   rb_max         upper bound on bi-elliptic apoapsis radius [km]
%   k_one_tangent  default one-tangent apsis factor (used when no TOF pressure)
%   weights        MCDA weights [dv, propellant, TOF, complexity]
m = struct('alt1', 400, 'alt2', 35786, 'apo1', NaN, 'apo2', NaN, 'dw_deg', 0, 'di_deg', 28.5, 'dv_budget', 6.5, 'margin_pct', 5, ...
           'tof_max_h', 24*200, 'propulsion', 'both', 'isp_chem', 320, 'tw_chem', 0.5, ...
           'isp_ep', 1600, 'accel_ep', 5e-4, 'rb_max', 384400, ...
           'k_one_tangent', 1.5, 'weights', [0.30 0.20 0.30 0.20], 'w1_deg', 0, 'elements', [], 'R', []);
for i = 1:2:numel(varargin)
    if ~isfield(m, varargin{i})
        error('odst_mission:badField', 'Unknown mission field "%s".', varargin{i});
    end
    m.(varargin{i}) = varargin{i+1};
end
if ~isempty(m.elements)
    [m.di_deg, m.w1_deg, m.dw_deg, m.R] = odst_relative_geometry(m);
end
end
