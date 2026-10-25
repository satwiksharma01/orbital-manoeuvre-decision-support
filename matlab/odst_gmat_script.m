function txt = odst_gmat_script(m, o, report)
%ODST_GMAT_SCRIPT NASA GMAT script that flies option o's impulsive burn plan.
%   txt = odst_gmat_script(m, o)          m: odst_mission, o: an option from odst_evaluate
%   txt = odst_gmat_script(m, o, report)  file the Arrival report is written to
%   The script starts from the tool's start state (o.y0), applies each burn of o.vplan as an
%   inertial impulse and reports the orbit it arrives on, with the tool's physics: point-mass
%   Earth with the tool's mu. A mission given by orbital elements is turned into
%   EarthMJ2000Eq (odst_inertial_plan). Same script as gmat_script in odst.py.
if nargin < 3, report = 'odst_arrival.txt'; end
C = odst_const();
a1 = m.alt1; if ~isnan(m.apo1), a1 = m.apo1; end
a2 = m.alt2; if ~isnan(m.apo2), a2 = m.apo2; end
[y0, vplan] = odst_inertial_plan(m, o);
n = size(vplan, 1);
if isempty(m.R)
    orbits = {sprintf('%% Start orbit %g x %g km altitude; target %g x %g km, di %g deg, dw %g deg', ...
                      m.alt1, a1, m.alt2, a2, m.di_deg, m.dw_deg)};
    frame = C.FRAME;
else
    e = m.elements;
    orbits = {sprintf('%% Start orbit %g x %g km altitude, i %g deg, RAAN %g deg, argp %g deg', m.alt1, a1, e(1:3))
              sprintf('%% Target orbit %g x %g km altitude, i %g deg, RAAN %g deg, argp %g deg; planes %.4g deg apart', ...
                      m.alt2, a2, e(4:6), m.di_deg)};
    frame = C.FRAME_ELEMENTS;
end
L = [{'% Orbital Manoeuvre Decision Support Tool: GMAT export'
      sprintf('%% Transfer: %s, %d burn(s), total dv %.4f km/s (impulsive), flight time %.2f h', ...
              o.name, n, o.dv, o.tof/3600)}; orbits
     {['% Frame: ' frame '.']
     '%   The epoch is arbitrary: the tool has no timing model yet, so Sat starts at the departure point.'
     '% Physics as in the tool: point-mass Earth with the tool''s mu and impulsive burns, so the'
     '%   Arrival report matches the target orbit to rounding. Edit TwoBody to add J2, drag or third bodies.'
     ''
     '% GMAT''s default is 398600.4415, enough to shift a 10-day transfer''s arrival by 0.2 km'
     ['Earth.Mu = ' num(C.MU) ';']
     ''
     'Create Spacecraft Sat;'
     'Sat.DateFormat = UTCGregorian;'
     ['Sat.Epoch = ''' C.EPOCH ''';']
     'Sat.CoordinateSystem = EarthMJ2000Eq;'
     'Sat.DisplayStateType = Cartesian;'}];
k = {'X', 'Y', 'Z', 'VX', 'VY', 'VZ'};
for j = 1:6, L{end+1} = ['Sat.' k{j} ' = ' num(y0(j)) ';']; end %#ok<AGROW>
L = [L; {''; 'Create ForceModel TwoBody;'; 'TwoBody.CentralBody = Earth;'; 'TwoBody.PointMasses = {Earth};'
         'Create Propagator Prop;'; 'Prop.FM = TwoBody;'; 'Prop.Type = RungeKutta89;'; 'Prop.InitialStepSize = 60;'
         'Prop.Accuracy = 1e-13;'; 'Prop.MinStep = 0;'; 'Prop.MaxStep = 600;'}];
for b = 1:n
    L = [L; {''; sprintf('Create ImpulsiveBurn Burn%d;', b); sprintf('Burn%d.CoordinateSystem = EarthMJ2000Eq;', b)}]; %#ok<AGROW>
    for j = 1:3, L{end+1} = sprintf('Burn%d.Element%d = %s;', b, j, num(vplan(b, j + 1))); end %#ok<AGROW>
end
L = [L; {''; 'Create ReportFile Arrival;'; ['Arrival.Filename = ''' report ''';']; 'Arrival.Precision = 16;'
         'Arrival.WriteHeaders = true;'; ''; 'Create OrbitView View;'; 'View.Add = {Sat, Earth};'; ''
         'BeginMissionSequence;'}];
t = 0;
for b = 1:n
    tb = vplan(b, 1);
    if tb > t, L{end+1} = ['Propagate Prop(Sat) {Sat.ElapsedSecs = ' num(tb - t) '};']; end %#ok<AGROW>
    L{end+1} = sprintf('Maneuver Burn%d(Sat);', b); %#ok<AGROW>
    t = tb;
end
r2 = C.RE + m.alt2; ra2 = C.RE + a2;
L{end+1} = ['Report Arrival ' strjoin(strcat('Sat.EarthMJ2000Eq.', k), ' ') ...
            ' Sat.Earth.SMA Sat.Earth.ECC Sat.EarthMJ2000Eq.INC Sat.EarthMJ2000Eq.RAAN Sat.EarthMJ2000Eq.AOP;'];
L{end+1} = ['Propagate Prop(Sat) {Sat.ElapsedSecs = ' num(2*pi*sqrt(((r2 + ra2)/2)^3/C.MU)) '};'];
txt = [strjoin(L, newline) newline];
end

function s = num(x)
% shortest text that reads back as the same double; + 0 turns -0 into 0
x = x + 0;
for p = 15:17
    s = sprintf('%.*g', p, x);
    if str2double(s) == x, return; end
end
end
