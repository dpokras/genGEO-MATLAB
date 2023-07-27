clear;
close all;
clc;
tic;

% First find the optimal m-dot. Then, use m to find Tprod
%% Define parameters for base case
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
%% IDEAL / BASELINE
params.wellCostType = 'Baseline';
%params.wellCostType = 'Ideal';
%% LW
params.res_length = 5000;
%params.res_length = 18435.279025; %0.25m radius, 3.5km, ideal
%params.res_length = 28986.7093; %0.25m radius, 8km, ideal
%params.res_length = 8487.45215; %9-5/8 inch diameter, 3.5km, ideal
%params.res_length = 12018.54335; %9-5/8 inch diameter, 8km, ideal
%params.res_length = 5895.0246; %0.25m radius, 3.5km, base case
%params.res_length = 8318.931705; %0.25m radius, 8km, base case
%params.res_length = 3899.643335; %9-5/8 inch diameter, 3.5km, base case
%params.res_length = 4930.954335; %9-5/8 inch diameter, 8km, base case

%% DIAMETER
%d = 0.244475; %9-5/8 inch in meters
d = 0.5; %base case
%d = 0.1016; % 4-inch in meters
params.well_radius = d/2;
%% DEPTH
params.depth = 3500;
%% ALL OTHER PARAMETERS
params.system = 'Conduction1';
params.wellFieldType = 'Doublet';
params.time_years = 30;
params.dT_dz = 0.035;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';
fluid_All = ["Water","CO2"];
% Select Fluid Type
params.fluid = fluid_All(1);
%% Calculate Outputs for Baseline
%params.m_dot_IP = 188; % random high number
%params.m_dot_IP = 22.6567; % CO2 opt m dot
params.m_dot_IP = 11.0297; % Water opt m dot
result = total_analytic_system_water(params);
T_out = result.T_prod_surface_C;
pump_pressure = result.dP_pump
W_pump = result.W_pump_inj_IP
P_out = result.productionUpper.EndPressure
Q_well = result.Q_fluid_IP
W_net = result.W_net_IP;
n_efficiency = W_net/Q_well;
result.injection.State(:,2)
result.reservoir.State(:,2)
result.productionLower.State(:,2)
result.productionUpper.State(:,2)
% output_baseline = [params.depth, params.res_length, result.m_dot_IP, result.W_net_IP/1e6, result.T_prod_surface_C, result.W_net_IP/(2*params.depth + 4*params.res_length), result.CapitalCost.C_greenfield,result.SpecificCapitalCost_greenfield];
% spcc = output_baseline(end)
% power = result.W_net_IP/1e6
%%% Calculate Outputs for Ideal
% params.wellCostType = 'Ideal';
% result = total_analytic_system_optmdot(params);
% output_ideal = [params.depth, params.res_length, result.m_dot_IP, result.W_net_IP/1e6, result.T_prod_surface_C, result.W_net_IP/(2*params.depth + 4*params.res_length), result.CapitalCost.C_greenfield, result.SpecificCapitalCost_greenfield]
% cc_prod = result.CapitalCost.C_wells_production/1e6;
% cc_inj = result.CapitalCost.C_wells_injection/1e6;
% cc_wellfield = result.CapitalCost.C_wellfield/1e6;
% cc_total = result.CapitalCost.C_greenfield/1e6
% percent = (cc_prod + cc_inj)*100 / cc_total;
% output = [output_baseline;output_ideal];
% %% Write Tables
% % For LCOE
% T1 = table(output);
% writetable(T1,'data\Base Case Results_water.xlsx')
