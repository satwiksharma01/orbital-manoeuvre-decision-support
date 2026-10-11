function res = odst_run_simulink(opt, r1, r2, diDeg, tau)
%ODST_RUN_SIMULINK Simulate a manoeuvre with finite burns in Simulink.
%   res = odst_run_simulink(opt, r1, r2, diDeg, tau)
%   opt   option struct from odst_evaluate (impulsive plan in opt.plan)
%   tau   burn duration per impulse [s] (default 60)
%   Compares the final orbit with the target (the impulsive plan reaches it exactly,
%   see odst_tests), so the residual is the finite-burn (gravity) loss.
if nargin < 5, tau = 60; end
if isempty(opt.plan)
    error('odst_run_simulink:noPlan', '%s has no impulsive burn plan to simulate.', opt.name);
end
C = odst_const(); mu = C.MU;
mdl = 'odst_transfer_sim';
if ~exist([mdl '.slx'], 'file'), odst_build_simulink(); end
load_system(mdl);

plan = opt.plan;
coast = 2*pi*sqrt(r2^3/mu);                       % one target-orbit period after arrival
bp = struct('t', plan(:, 1), 'vafter', plan(:, 2), 'dinc', plan(:, 3), 'tau', tau, 'mu', mu);
n1 = sqrt(mu/r1^3); th0 = -n1*tau;                % start tau earlier on the initial orbit
v1 = sqrt(mu/r1);
x0 = [r1*cos(th0); r1*sin(th0); 0; -v1*sin(th0); v1*cos(th0); 0];
assignin('base', 'odst_burnplan', bp);
assignin('base', 'odst_x0', x0);
assignin('base', 'odst_tstop', plan(end, 1) + tau + coast);
clear odst_sim_dynamics

out = sim(mdl);
T = out.get('simT'); X = out.get('simX');
[a, e, inc] = odst_elements(X(end, :).');
res = struct('a', a, 'e', e, 'inc', inc, 'a_err', a - r2, 'inc_err', inc - diDeg);
fprintf('\nSimulink finite-burn verification: %s (tau = %g s per burn)\n', opt.name, tau);
fprintf('  final a   = %.2f km  (target %.2f, error %+.2f km)\n', a, r2, a - r2);
fprintf('  final e   = %.2e\n', e);
fprintf('  final inc = %.3f deg (target %.3f, error %+.3f deg)\n', inc, diDeg, inc - diDeg);

TH = odst_theme();                          % shared design tokens (T is the time vector here)
lineColor = TH.series.(opt.key);
fig = figure('Name', 'Simulink verification', 'Position', [100 100 1100 430]);
if isprop(fig, 'Theme'), fig.Theme = 'light'; end
fig.Color = 'w';
subplot(1, 2, 1); hold on; axis equal; grid on;
th = linspace(0, 2*pi, 400);
k = 1e-3;                                   % plot in 10^3 km (no axis exponent)
fill(k*C.RE*cos(th), k*C.RE*sin(th), TH.earth, 'EdgeColor', 'none');
plot(k*r1*cos(th), k*r1*sin(th), 'Color', TH.orbit, 'LineWidth', 1);
plot(k*r2*cos(th), k*r2*sin(th), '--', 'Color', TH.orbit, 'LineWidth', 1);
plot(k*X(:, 1), k*X(:, 2), 'LineWidth', 2, 'Color', lineColor);
xlabel('x [10^3 km]'); ylabel('y [10^3 km]');
title(['Simulink: ' opt.name]);
subtitle(sprintf('Finite burns of \\tau = %g s replace each impulse', tau));
subplot(1, 2, 2); hold on; grid on;
rn = sqrt(sum(X(:, 1:3).^2, 2));
plot(T/3600, k*rn, 'LineWidth', 2, 'Color', lineColor);
yline(k*r2, '--', 'target r_2');
xlabel('time [h]'); ylabel('|r| [10^3 km]');
title('Radius history');
subtitle(sprintf('Final a error %+.2f km, e = %.1e, i error %+.3f\\circ', a - r2, e, inc - diDeg));
end
