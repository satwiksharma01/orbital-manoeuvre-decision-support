function [names, W] = odst_presets()
%ODST_PRESETS Mission-profile weight presets [dv, propellant, TOF, complexity].
names = {'Balanced', 'Propellant-critical', 'Time-critical', 'Low-complexity'};
W = [0.30 0.20 0.30 0.20;    % Balanced
     0.25 0.55 0.05 0.15;    % Propellant-critical (e.g. GEO comsat, all-electric raising)
     0.15 0.10 0.60 0.15;    % Time-critical (crew, rapid response, debris avoidance)
     0.20 0.20 0.10 0.50];   % Low-complexity (CubeSat / limited autonomy)
end
