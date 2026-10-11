function mdl = odst_build_simulink()
%ODST_BUILD_SIMULINK Programmatically build the ODST finite-burn Simulink model.
%   Creates and saves odst_transfer_sim.slx next to this file:
%
%     Clock --+------------------------------> To Workspace (simT)
%             |
%             v
%     Mux [t; x] --> Interpreted MATLAB Fcn --> Integrator (x = [r; v]) --+--> To Workspace (simX)
%       ^             odst_sim_dynamics(u)       IC = odst_x0             |
%       +-----------------------------------------------------------------+
%
%   odst_sim_dynamics returns xdot = [v; -mu r/|r|^3 + a_thrust(t)], where the
%   thrust acceleration reproduces each impulsive burn of the selected plan as a
%   constant-acceleration finite burn of duration tau centred on the burn time.
%   Run it with odst_run_simulink.
mdl = 'odst_transfer_sim';
folder = fileparts(mfilename('fullpath'));
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
open_system(mdl);

add_block('simulink/Sources/Clock', [mdl '/Clock'], 'Position', [40 60 70 90]);
add_block('simulink/Signal Routing/Mux', [mdl '/Mux'], 'Inputs', '2', ...
          'Position', [140 60 145 140]);
add_block('simulink/User-Defined Functions/Interpreted MATLAB Function', [mdl '/Dynamics'], ...
          'MATLABFcn', 'odst_sim_dynamics(u)', 'OutputDimensions', '6', ...
          'Position', [200 80 330 120]);
add_block('simulink/Continuous/Integrator', [mdl '/State'], ...
          'InitialCondition', 'odst_x0', 'Position', [390 85 420 115]);
add_block('simulink/Sinks/To Workspace', [mdl '/ToX'], 'VariableName', 'simX', ...
          'SaveFormat', 'Array', 'Position', [500 85 560 115]);
add_block('simulink/Sinks/To Workspace', [mdl '/ToT'], 'VariableName', 'simT', ...
          'SaveFormat', 'Array', 'Position', [500 20 560 50]);
add_block('simulink/Sinks/Scope', [mdl '/Scope'], 'Position', [500 150 530 180]);

add_line(mdl, 'Clock/1', 'Mux/1', 'autorouting', 'on');
add_line(mdl, 'State/1', 'Mux/2', 'autorouting', 'on');
add_line(mdl, 'Mux/1', 'Dynamics/1', 'autorouting', 'on');
add_line(mdl, 'Dynamics/1', 'State/1', 'autorouting', 'on');
add_line(mdl, 'State/1', 'ToX/1', 'autorouting', 'on');
add_line(mdl, 'State/1', 'Scope/1', 'autorouting', 'on');
add_line(mdl, 'Clock/1', 'ToT/1', 'autorouting', 'on');

set_param(mdl, 'Solver', 'ode45', 'RelTol', '1e-9', 'AbsTol', '1e-9', ...
          'StartTime', '-odst_burnplan.tau', 'StopTime', 'odst_tstop', ...
          'MaxStep', 'odst_burnplan.tau/4', 'ReturnWorkspaceOutputs', 'on');
save_system(mdl, fullfile(folder, [mdl '.slx']));
fprintf('Saved %s\n', fullfile(folder, [mdl '.slx']));
end
