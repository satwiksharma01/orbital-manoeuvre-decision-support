function C = odst_const()
%ODST_CONST Physical constants used by the Orbital Manoeuvre Decision Support Tool.
%   Units: km, s, km/s.
C.MU = 398600.4418;   % Earth gravitational parameter [km^3/s^2]
C.RE = 6378.137;      % Earth equatorial radius [km]
C.G0 = 9.80665e-3;    % standard gravity [km/s^2]
% GMAT and JSON export (odst_gmat_script, odst_burns_json): the epoch is arbitrary, as
% the tool has no timing model yet
C.EPOCH = '01 Jan 2026 12:00:00.000';    % UTC
C.FRAME = 'EarthMJ2000Eq, start orbit in the X-Y plane with its periapsis on +X';
% Operational complexity index per manoeuvre (1 = simplest ... 5 = hardest).
% See docs/REPORT.md section 4.4 for the justification of each value.
C.CPLX = struct('hohmann', 1.0, 'hohmann_split', 2.0, 'one_tangent', 2.5, ...
                'bielliptic', 3.0, 'low_thrust', 4.0, 'apse', 1.5, 'lambert', 2.5);
end
