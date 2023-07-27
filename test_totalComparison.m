clear all;
close all;

params = SimulationParameters;
params.system = 'Porous';
params.optimizationMode = 'MinimizeLCOE_Brownfield';
params.orcFluid = 'R245fa';
params.res_length = 707;
params.F_OM = 0.055;
params.discountRate = 0.04;        
params.CapacityFactor = 0.95;
params.time_years = 1;
params.modelResTemperatureDepletion = false;
params.well_radius = 0.155;



%2.5km, CO2
% params.fluid = 'CO2';
% params.depth = 2500; 
% params.dT_dz = 0.035;
% params.transmissivity =  1e-15 * 50 * 300;
% params.wellFieldType = '5spot_SharedNeighbor';
% params.simulatorType = 'genGEO';
% 
% result_baseCase_co2 = total_analytic_system_optmdot(params);
% 
% %2.5km, Water
% params.fluid = 'Water';
% 
% result_baseCase_water = total_analytic_system_optmdot(params);
% 
% %2.5km, Geophires
% params.simulatorType = 'GEOPHIRES';
% result_baseCase_geophires = total_analytic_system_optmdot(params);
% 
% 
% 
% 
% %4km, CO2
% params.fluid = 'CO2';
% params.depth = 3630; 
% params.dT_dz = 0.035;
% params.transmissivity =  1e-15 * 50 * 300;
% params.wellFieldType = '5spot_SharedNeighbor';
% params.simulatorType = 'genGEO';
% 
% result_4km_co2 = total_analytic_system_optmdot(params);
% 
% %4km, Water
% params.fluid = 'Water';
% result_4km_water = total_analytic_system_optmdot(params);
% 
% %4km, Geophires
% params.simulatorType = 'GEOPHIRES';
% result_4km_geophires = total_analytic_system_optmdot(params);



%Tungsten, CO2
params.fluid = 'CO2';
params.depth = 1000; 
params.dT_dz = 0.127;
params.transmissivity =  1e-9;
params.wellFieldType = 'Tungsten';
params.simulatorType = 'genGEO';

result_tungsten_co2 = total_analytic_system_optmdot(params);

%Tungsten, Water
params.m_dot_IP = 250;
params.fluid = 'Water';
result_tungsten_water = total_analytic_system(params);

% Tungsten, Geophires
params.simulatorType = 'GEOPHIRES';
result_tungsten_geophires = total_analytic_system(params);

%
Lastly, Tungs



