clear;
clc;

params                  = SimulationParameters; %radius, surface temp, dT/dz, etc.

params.time_years       = 1:30;
params.wellCostType     = 'Baseline';
% params.wellCostType     = 'Ideal';
params.find_opt_S_ratio_toggle = 1;
params.S_ratio = 0.5; %S_ratio = m_S1/m_total
params.config = 1;
if params.config == 1 || params.config == 2
    params.S_ratio = 1;
end

params.well_radius = 0.4667/2; %m ~18 3/8" inner diameter
params.side_stream_radius = 0.2032/2; %m ~8" inner diameter


params.coolingMode = 'Wet';
params.n_streams = 10;
depths = params.depth - params.thickness*((params.n_streams-1):-1:0);
temp_at_depth = params.T_surface_rock_C + params.dT_dz * depths; %C
tic;
result_injWell = semi_analytic_well(1e7, 30, params.T_surface_rock_C, -depths, 0, 70, params);
toc;
a=1