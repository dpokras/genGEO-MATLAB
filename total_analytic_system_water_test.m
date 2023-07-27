clear all;
close all;

params = SimulationParameters;
params.depth = 2500;
params.transmissivity = 1.0e-15 * 15000;
params.thickness = 100;
params.time_years = 1;
params.fluid = 'Water';
params.system = 'Porous';
params.CapacityFactor = 0.9;


%Test1 Small m_dot = 1
params.m_dot_IP = 1;
result = total_analytic_system_water(params);

Assert(result.T_prod_surface_C, 81.2595, 'test1_T_prod_surface_C');
Assert(result.W_net_IP, 5.2775e3, 'test1_W_net_IP');
Assert(result.CapitalCost.C_brownfield, 9.1965e6, 'test1_C_brownfield_N');
Assert(result.CapitalCost.C_greenfield, 2.2308e7, 'test1_C_greenfield_N');
Assert(result.LCOE_brownfield, 0.0083879, 'test1_LCOE_brownfield');

%Test2 Medium m_dot = 40
params.m_dot_IP = 40;
result = total_analytic_system_water(params);

Assert(result.T_prod_surface_C, 100.3125, 'test2_T_prod_surface_C');
Assert(result.W_net_IP, 1.4286e+05, 'test2_W_net_IP');
Assert(result.CapitalCost.C_brownfield, 2.3287e7, 'test2_C_brownfield_N');
Assert(result.CapitalCost.C_greenfield, 3.6399062e7, 'test2_C_greenfield_N');
Assert(result.LCOE_brownfield, 7.8465e-4, 'test2_LCOE_brownfield');

    