clear;
close all;clc
% Purpose of sheet is to generate a plot for the different system
% configurations (1 Lateral, 2 Laterals,...) for a variety of reservoir
% depths
%% Define System Parameters
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
params.wellFieldType = 'Doublet'; 
params.time_years = 30;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';
params.depth = 8000;

% Select Fluid Type
fluid_All = ["Water","CO2"];
params.fluid = fluid_All(2);
%% Base Case
params.dT_dz = 0.035;
params.system = "Conduction4";
params.res_length = 5000;
params.well_radius = 0.25;
params.wellCostType = 'Baseline';
params.wellCostType = 'Ideal';

%% Num_Laterals
% params.system = "Conduction8";
% params.res_length = 4174.1908965;
% 
%% Diameter
% params.well_radius = 0.244475/2;
% params.res_length = 3899.643335;

%% Gradient
params.dT_dz = 0.040;
params.res_length = 5997.680148;

result = total_analytic_system_optmdot(params);
Power = result.W_net_IP / 1e3 % in kW
LCOE = result.LCOE_greenfield * 1e6 %in 2019$/MWh