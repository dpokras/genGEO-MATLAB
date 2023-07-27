clear all;
close all;

params = SimulationParameters;
params.depth = 2500;
params.transmissivity = 1.0e-15 * 15000;
params.thickness = 100;
params.time_years = 1;
params.fluid = 'CO2';
params.system = 'Porous';
params.optimizationMode = 'MinimizeLCOE_Brownfield';
params.CapacityFactor = 0.9;



%Test1 Small m_dot = 10
params.m_dot_IP = 10;
result = total_analytic_system_co2(params);

Assert(result.dP_turbine, 5.9421e6, 'test1_dP_turbine');
Assert(result.T_prod_surface_C, 59.2802, 'test1_T_prod_surface_C');
Assert(result.W_net_IP, 9.7662e4, 'test1_W_net_IP');
Assert(result.CapitalCost.C_brownfield, 1.3271e7, 'test1_C_brownfield_N');
Assert(result.CapitalCost.C_greenfield, 3.5901e7, 'test1_C_greenfield_N');
Assert(result.LCOE_brownfield, 6.5410e-4, 'test1_LCOE_brownfield');

%Test2 Medium m_dot = 80
params.m_dot_IP = 80;
result = total_analytic_system_co2(params);

Assert(result.dP_turbine, 5.6501e6, 'test2_dP_turbine');
Assert(result.T_prod_surface_C, 59.0540, 'test2_T_prod_surface_C');
Assert(result.W_net_IP, 4.9940e+05, 'test2_W_net_IP');
Assert(result.CapitalCost.C_brownfield, 2.7650e7, 'test2_C_brownfield_N');
Assert(result.CapitalCost.C_greenfield, 5.0280e7, 'test2_C_greenfield_N');
Assert(result.LCOE_brownfield, 2.6650e-4, 'test2_LCOE_brownfield');
    

%Test3 Large m_dot = 200
params.m_dot_IP = 200;
result = total_analytic_system_co2(params);

Assert(result.dP_turbine, 3.468765e+06, 'test3_dP_turbine');
Assert(result.T_prod_surface_C, 47.5801, 'test3_T_prod_surface_C');
Assert(result.W_net_IP, -1.3602e6, 'test3_W_net_IP');
Assert(result.CapitalCost.C_brownfield, 5.0443e7, 'test3_C_brownfield_N');
Assert(result.CapitalCost.C_greenfield, 7.3073e7, 'test3_C_greenfield_N');
Assert(result.LCOE_brownfield, NaN, 'test3_LCOE_brownfield');

%Test4 Throttle, not pump
params.m_dot_IP = 100;
params.depth = 2400;
params.transmissivity = 1e-8;
result = total_analytic_system_co2(params);

Assert(result.dP_turbine, 5.0019e+06, 'test4_dP_turbine');
Assert(result.dP_pump, -1.0756e6, 'test4_dP_pump');
Assert(result.T_prod_surface_C, 55.3144, 'test4_T_prod_surface_C');
Assert(result.W_net_IP, 8.2656e5, 'test4_W_net_IP');
Assert(result.LCOE_brownfield, 1.5247e-4, 'test4_LCOE_brownfield');
    