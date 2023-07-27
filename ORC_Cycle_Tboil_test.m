clear all;
close all;

params = SimulationParameters;
params.coolingMode = 'Wet';
params.orcFluid = 'R245fa';

T_in_C = 150;
T_boil_C = 100;

%Test1 Regular
result = ORC_Cycle_Tboil(T_in_C, T_boil_C, params.dT_orc_pinch, params);
    
Assert(result.T_out_C, 68.36, 'test1_temp');
Assert(result.w_net, 3.8559e4, 'test1_w_net');
Assert(result.w_turbine, 4.7773e4, 'test1_w_turbine');
Assert(result.q_preheater, 1.5778e5, 'test1_q_preheater');
Assert(result.q_boiler, 1.9380e5, 'test1_q_boiler');

%Test2, try with lower T_in than boiling
T_in_C = 15;

try
    result = ORC_Cycle_Tboil(T_in_C, T_boil_C, params.dT_orc_pinch, params);
    success = true;
catch
    success = false;
end
    
Assert(success, false, 'test2');
