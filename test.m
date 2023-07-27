clear all;
close all;

tic;

params = SimulationParameters;
params.m_dot_IP = 100;
params.depth = 2500; 
%params.dT_dz = 0.127;
%params.well_radius = 0.122;
params.res_length = 707;
params.transmissivity =  1e-15 * 50 * 300;
%params.transmissivity =  1e-9;
params.time_years = 1;
%params.time_years = 30;
params.fluid = 'CO2';
%params.fluid = 'Water';
params.system = 'Porous';
%params.system = 'Conduction1';
%params.system = 'Conduction2';
%params.system = 'Conduction4';
%params.system = 'Conduction8';
params.optimizationMode = 'MinimizeLCOE_Brownfield';
%params.optimizationMode = 'MinimizeLCOE_Greenfield';
%params.optimizationMode = 'MaximizePower';
params.orcFluid = 'R245fa';
%params.orcFluid = 'R600a';
%params.wellFieldType = 'Doublet';
%params.wellFieldType = 'Tungsten';
params.wellFieldType = '5spot_SharedNeighbor';
%params.wellCostType = 'Ideal';
%params.hasSurfaceGatheringSystem = false;
params.modelResTemperatureDepletion = false;
params.F_OM = 0.055;
params.discountRate = 0.04;        
params.CapacityFactor = 0.95;
%params.coolingMode = 'Dry';
%params.dT_approach = 15;

result = total_analytic_system(params);
%result = total_analytic_system_optmdot(params);



W_net_IP_MW = result.W_net_IP / 1e6;
Q_fluid_IP_MW = result.Q_fluid_IP / 1e6;
dP_turbine_MPa = result.dP_turbine / 1e6;
SpecificCapitalCost_greenfield = result.SpecificCapitalCost_greenfield;
LCOE_greenfield_MWh = result.LCOE_greenfield * 1e6; 
SpecificCapitalCost_brownfield = result.SpecificCapitalCost_brownfield;
LCOE_brownfield_MWh = result.LCOE_brownfield * 1e6;



toc;