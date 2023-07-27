clear all;
close all;

Psi_1 = 2;
p1 = -1.9376;
p2 = 1.743;
p3 = 0.182;
[gamma] = DepletionCurve(Psi_1, p1, p2, p3);

Assert(gamma, 0.2531, 'test1');


params = SimulationParameters;
params.thickness = 300;
permeability = 50e-15;
params.transmissivity = params.thickness * permeability;
params.depth = 2500;
params.dT_dz = 0.035;
params.well_radius = 0.205;
params.time_years = 30;
params.fluid = 'CO2';

P_f_initial = 25e6;
T_f_initial = 40;
m_dot = 100;
params.wellFieldType = '5spot';

%Test2, no transients
params.modelResPressureTransient = false;
params.modelResTemperatureDepletion = false;

result = PorousReservoir(P_f_initial, T_f_initial, m_dot, params);

Assert(result.EndPressure, 2.0351e7, 'test2_pressure');
Assert(result.EndTemp, 102.5, 'test2_temp');
Assert(result.EndEnthalpy, 4.30e5, 'test2_enthalpy');
Assert(result.Psi, 1.0151, 'test2_Psi');

%Test3, temp depletion
params.modelResPressureTransient = false;
params.modelResTemperatureDepletion = true;

result = PorousReservoir(P_f_initial, T_f_initial, m_dot, params);

Assert(result.EndPressure, 2.0351e7, 'test3_pressure');
Assert(result.EndTemp, 71.7201, 'test3_temp');
Assert(result.EndEnthalpy, 3.5324e5, 'test3_enthalpy');
Assert(result.Psi, 0.8457, 'test3_Psi');

%Test4, temp + pressure change
params.modelResPressureTransient = true;
params.modelResTemperatureDepletion = true;

result = PorousReservoir(P_f_initial, T_f_initial, m_dot, params);

Assert(result.EndPressure, 1.932e7, 'test4_pressure');
Assert(result.EndTemp, 71.7201, 'test4_temp');
Assert(result.EndEnthalpy, 3.5732e5, 'test4_enthalpy');
Assert(result.Psi, 0.8457, 'test4_Psi');

%Test5, temp depletion, higher flowrate
m_dot = 200;
params.modelResPressureTransient = false;
params.modelResTemperatureDepletion = true;

result = PorousReservoir(P_f_initial, T_f_initial, m_dot, params);

Assert(result.EndPressure, 1.5702e7, 'test5_pressure');
Assert(result.EndTemp, 55.8873, 'test5_temp');
Assert(result.EndEnthalpy, 3.2876e5, 'test5_enthalpy');
Assert(result.Psi, 1.2962, 'test5_Psi');

