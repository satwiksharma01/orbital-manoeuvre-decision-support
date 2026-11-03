function txt = odst_gmat_script(m, o, report)
%ODST_GMAT_SCRIPT NASA GMAT script that flies option o's impulsive burn plan.
%   txt = odst_gmat_script(m, o)          m: odst_mission, o: an option from odst_evaluate
%   txt = odst_gmat_script(m, o, report)  file the Arrival report is written to
%   The script starts from the tool's start state (o.y0), applies each burn of o.vplan as an
%   inertial impulse and reports the orbit it arrives on, with the tool's physics: point-mass
%   Earth with the tool's mu. A mission given by orbital elements is turned into
%   EarthMJ2000Eq (odst_inertial_plan). With J2 (m.j2): a J2 field (JGM-2), each burn in its
%   local VNB frame at its apsis, time or node (odst_local_burns), and a TargetSat flying the
%   target orbit alongside. Same script as gmat_script in odst.py.
if nargin < 3, report = 'odst_arrival.txt'; end
C = odst_const();
a1 = m.alt1; if ~isnan(m.apo1), a1 = m.apo1; end
a2 = m.alt2; if ~isnan(m.apo2), a2 = m.apo2; end
j2 = m.j2 && ~isempty(m.elements);
if j2
    [y0, burns] = odst_local_burns(m, o);
else
    [y0, vplan] = odst_inertial_plan(m, o);
    burns = struct('t', num2cell(vplan(:, 1)), 'dv', num2cell(vplan(:, 2:4), 2), 'at', 'time', 'signed', false);
end
n = numel(burns);
if j2
    physics = {'% Physics as in the tool''s J2 plan: Earth''s J2 (JGM-2, degree 2, order 0) with the tool''s mu. Each burn'
               '%   is in its local VNB frame (velocity, orbit normal, binormal), so it follows the plane J2 turns, and fires'
               '%   at its apsis, time or equator crossing. TargetSat flies the target orbit alongside: compare their planes.'
               '%   J2''s effects within an orbit, which the tool leaves out, leave a small mismatch for a trim burn.'};
else
    physics = {'% Physics as in the tool: point-mass Earth with the tool''s mu and impulsive burns, so the'
               '%   Arrival report matches the target orbit to rounding. Edit TwoBody to add J2, drag or third bodies.'};
end
wait = '';
if strcmp(o.key, 'nodal_drift')
    wait = sprintf(', including %.2f days at %.1f km', o.params.wait_days, o.params.alt_d);
end
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
      sprintf('%% Transfer: %s, %d burn(s), total dv %.4f km/s (impulsive), flight time %.2f h%s', ...
              o.name, n, o.dv, o.tof/3600, wait)}; orbits
     {['% Frame: ' frame '.']
     '%   The epoch is arbitrary: the tool has no timing model yet, so Sat starts at the departure point.'}; physics
     {''
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
r2 = C.RE + m.alt2; ra2 = C.RE + a2;
if j2
    e = m.elements;
    L = [L; {''; 'Create Spacecraft TargetSat;'; 'TargetSat.DateFormat = UTCGregorian;'; ['TargetSat.Epoch = ''' C.EPOCH ''';']
             'TargetSat.CoordinateSystem = EarthMJ2000Eq;'; 'TargetSat.DisplayStateType = Keplerian;'
             ['TargetSat.SMA = ' num((r2 + ra2)/2) ';']; ['TargetSat.ECC = ' num((ra2 - r2)/(ra2 + r2)) ';']
             ['TargetSat.INC = ' num(e(4)) ';']; ['TargetSat.RAAN = ' num(mod(e(5), 360)) ';']
             ['TargetSat.AOP = ' num(mod(e(6), 360)) ';']; 'TargetSat.TA = 0;'}];
end
sc = ~isnan(m.mass0_kg);
if sc                     % GMAT draws each burn's propellant from the tank: Sat.Fuel.FuelMass is reported
    L = [L; {''; '% The spacecraft: the electric propellant is counted in the dry mass, as only chemical burns fly here'
             'Create ChemicalTank Fuel;'; ['Fuel.FuelMass = ' num(m.prop_chem_kg) ';']; 'Fuel.AllowNegativeFuelMass = true;'
             'Fuel.FuelDensity = 1000;'; ['Fuel.Volume = ' num(m.prop_chem_kg/1000) ';']
             ['Sat.DryMass = ' num(m.mass0_kg - m.prop_chem_kg) ';']; 'Sat.Tanks = {Fuel};'}];
end
if j2
    fm = 'EarthJ2'; sats = 'Sat, TargetSat';
    grav = {[fm '.PrimaryBodies = {Earth};']; [fm '.GravityField.Earth.Degree = 2;']; [fm '.GravityField.Earth.Order = 0;']
            [fm '.GravityField.Earth.PotentialFile = ''JGM2.cof'';']};
else
    fm = 'TwoBody'; sats = 'Sat'; grav = {[fm '.PointMasses = {Earth};']};
end
L = [L; {''; ['Create ForceModel ' fm ';']; [fm '.CentralBody = Earth;']}; grav
     {'Create Propagator Prop;'; ['Prop.FM = ' fm ';']; 'Prop.Type = RungeKutta89;'; 'Prop.InitialStepSize = 60;'
      'Prop.Accuracy = 1e-13;'; 'Prop.MinStep = 0;'; 'Prop.MaxStep = 600;'}];
for b = 1:n
    names = {sprintf('Burn%d', b)}; dvs = {burns(b).dv};
    if burns(b).signed, names{2} = sprintf('Burn%dD', b); dvs{2} = burns(b).dv.*[1, -1, 1]; end
    for q = 1:numel(names)
        nm = names{q};
        if j2
            L = [L; {''; ['Create ImpulsiveBurn ' nm ';']; [nm '.CoordinateSystem = Local;']; [nm '.Origin = Earth;']; [nm '.Axes = VNB;']}]; %#ok<AGROW>
        else
            L = [L; {''; ['Create ImpulsiveBurn ' nm ';']; [nm '.CoordinateSystem = EarthMJ2000Eq;']}]; %#ok<AGROW>
        end
        for j = 1:3, L{end+1} = sprintf('%s.Element%d = %s;', nm, j, num(dvs{q}(j))); end %#ok<AGROW>
        if sc
            L = [L; {[nm '.DecrementMass = true;']; [nm '.Tank = {Fuel};']
                     [nm '.Isp = ' num(m.isp_chem) ';']; [nm '.GravitationalAccel = ' num(C.G0*1e3) ';']}]; %#ok<AGROW>
        end
    end
end
L = [L; {''; 'Create ReportFile Arrival;'; ['Arrival.Filename = ''' report ''';']; 'Arrival.Precision = 16;'
         'Arrival.WriteHeaders = true;'; ''; 'Create OrbitView View;'; ['View.Add = {' sats ', Earth};']; ''
         'BeginMissionSequence;'}];
t = 0;
for b = 1:n
    tb = burns(b).t; at = burns(b).at;
    if any(strcmp(at, {'apoapsis', 'periapsis'}))
        L{end+1} = sprintf('Propagate Prop(%s) {Sat.Earth.%s};', sats, [upper(at(1)) at(2:end)]); %#ok<AGROW>
    elseif tb > t
        L{end+1} = ['Propagate Prop(' sats ') {Sat.ElapsedSecs = ' num(tb - t) '};']; %#ok<AGROW>
    end
    if strcmp(at, 'node'), L{end+1} = ['Propagate Prop(' sats ') {Sat.EarthMJ2000Eq.Z = 0};']; end %#ok<AGROW>
    if burns(b).signed
        L = [L; {'If Sat.EarthMJ2000Eq.VZ > 0'; sprintf('   Maneuver Burn%d(Sat);', b); 'Else'
                 sprintf('   Maneuver Burn%dD(Sat);', b); 'EndIf;'}]; %#ok<AGROW>
    else
        L{end+1} = sprintf('Maneuver Burn%d(Sat);', b); %#ok<AGROW>
    end
    t = tb;
end
L{end+1} = ['Report Arrival ' strjoin(strcat('Sat.EarthMJ2000Eq.', k), ' ') ...
            ' Sat.Earth.SMA Sat.Earth.ECC Sat.EarthMJ2000Eq.INC Sat.EarthMJ2000Eq.RAAN Sat.EarthMJ2000Eq.AOP' ...
            repmat(' Sat.Fuel.FuelMass', 1, sc) ...
            repmat(' TargetSat.EarthMJ2000Eq.INC TargetSat.EarthMJ2000Eq.RAAN TargetSat.Earth.SMA TargetSat.Earth.ECC TargetSat.EarthMJ2000Eq.AOP', 1, j2) ';'];
L{end+1} = ['Propagate Prop(' sats ') {Sat.ElapsedSecs = ' num(2*pi*sqrt(((r2 + ra2)/2)^3/C.MU)) '};'];
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
