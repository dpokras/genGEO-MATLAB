clear all;
close all;

params = SimulationParameters;
params.coolingMode = 'Wet';
params.orcFluid = 'R245fa';


%Test1 Regular
T_in_C = 150;
params.optimizationMode = 'MinimizeLCOE_Brownfield';
result = ORC_Cycle(T_in_C, params);
    
Assert(result.T_boil_C, 98.6284, 'test1_boil_temp');
Assert(result.T_out_C, 73.2095, 'test1_out_temp');
Assert(result.w_net, 3.5806e4, 'test1_w_net');
Assert(result.w_turbine, 4.4450e4, 'test1_w_turbine');
Assert(result.q_preheater, 1.4600e5, 'test1_q_preheater');
Assert(result.q_boiler, 1.8472e5, 'test1_q_boiler');

%Test2, try with lower T_in than boiling
T_in_C = 40;
params.optimizationMode = 'MinimizeLCOE_Brownfield';
result = ORC_Cycle(T_in_C, params);
    
Assert(result.T_boil_C, 30.1205, 'test2_boil_temp');
Assert(result.T_out_C, 34.8435, 'test2_out_temp');
Assert(result.w_net, -116.7655, 'test2_w_net');
Assert(result.w_turbine, 452.2658, 'test2_w_turbine');
Assert(result.q_preheater, 1.1579e3, 'test2_q_preheater');
Assert(result.q_boiler, 2.0393e4, 'test2_q_boiler');


%Test3 100C, max power
T_in_C = 100;
params.optimizationMode = 'MaximizePower';
result = ORC_Cycle(T_in_C, params);
    
Assert(result.T_boil_C, 64.3975, 'test3_boil_temp');
Assert(result.T_out_C, 58.7690, 'test3_out_temp');
Assert(result.w_net, 1.1073e4, 'test3_w_net');
Assert(result.w_turbine, 1.5493e4, 'test3_w_turbine');
Assert(result.q_preheater, 4.4804e4, 'test3_q_preheater');
Assert(result.q_boiler, 1.2900e5, 'test3_q_boiler');

%Test R600, max power
%Test4 Regular
T_in_C = 150;
params.orcFluid = 'R600a';
params.optimizationMode = 'MaximizePower';
result = ORC_Cycle(T_in_C, params);
    
Assert(result.T_boil_C, 99.2618, 'test4_boil_temp');
Assert(result.T_out_C, 59.4386, 'test4_out_temp');
Assert(result.w_net, 4.1002e4, 'test4_w_net');
Assert(result.w_turbine, 5.2823e4, 'test4_w_turbine');
Assert(result.q_preheater, 1.9304e5, 'test4_q_preheater');
Assert(result.q_boiler, 1.9698e5, 'test4_q_boiler');

%Test5 100C, min cost
T_in_C = 100;
params.optimizationMode = 'MinimizeLCOE_Brownfield';
result = ORC_Cycle(T_in_C, params);
    
Assert(result.T_boil_C, 66.5700, 'test5_boil_temp');
Assert(result.T_out_C, 61.4876, 'test5_out_temp');
Assert(result.w_net, 1.0661e4, 'test5_w_net');
Assert(result.w_turbine, 1.5112e4, 'test5_w_turbine');
Assert(result.q_preheater, 4.7441e4, 'test5_q_preheater');
Assert(result.q_boiler, 1.1491e5, 'test5_q_boiler');





