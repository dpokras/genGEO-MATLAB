clear;
close all;

params = SimulationParameters;
params.depth = 2500;
params.transmissivity = 1.0e-15 * 15000;
params.thickness = 100;
params.time_years = 1;
params.system = 'Porous';
params.optimizationMode = 'MinimizeLCOE_Brownfield';
params.CapacityFactor = 0.9;


%Test1 CO2
params.fluid = 'CO2';
tic
result = total_analytic_system_optmdot(params);
toc
% ~26 seconds to complete

Assert(result.m_dot_IP, 54.8493, 'test1_m_dot_IP');
Assert(result.W_net_IP, 4.5411e5, 'test1_W_net_IP');
Assert(result.LCOE_brownfield, 2.3915e-4, 'test1_LCOE_brownfield');


%Test2 Water
params.fluid = 'Water';
tic
result = total_analytic_system_optmdot(params);
toc
% ~35 seconds to complete

Assert(result.m_dot_IP, 21.9823, 'test2_m_dot_IP');
Assert(result.W_net_IP, 1.5212e5, 'test2_W_net_IP');
Assert(result.LCOE_brownfield, 5.4633e-4, 'test2_LCOE_brownfield');






    