clear all;
close all;

params = SimulationParameters;
params.dT_dz = 0.06;
params.time_years = 10;
params.fluid = 'Water';

P_f_initial = 1e6;
T_f_initial = 25;
dz_total_h = 0;
dr_total_h = 3000;
dz_total_v = -3500;
dr_total_v = 0;
m_dot = 5;
wellRadius = 0.279;

% Vertical well segment
T_e_initial = 15;
result = semi_analytic_well(P_f_initial, T_f_initial, T_e_initial, dz_total_v, dr_total_v, m_dot, wellRadius, params);

Assert(result.EndPressure, 3.533e7, 'VerticalInjection_Water_Pressure');
Assert(result.EndTemp, 67.03, 'VerticalInjection_Water_Temp');
Assert(result.EndEnthalpy, 3.0963e5, 'VerticalInjection_Water_Enthalpy');

% Horizontal well segment
T_res = T_e_initial + (params.dT_dz * abs(dz_total_v));
result = semi_analytic_well(result.EndPressure, result.EndTemp, T_res, dz_total_h, dr_total_h, m_dot, wellRadius, params);

Assert(result.EndPressure, 3.533e7, 'HorizontalInjection_Water_Pressure');
Assert(result.EndTemp, 121.99, 'HorizontalInjection_Water_Temp');
Assert(result.EndEnthalpy, 5.3712e5, 'HorizontalInjection_Water_Enthalpy');

% Now CO2
params.fluid = 'CO2';
% Vertical well segment
result = semi_analytic_well(P_f_initial, T_f_initial, T_e_initial, dz_total_v, dr_total_v, m_dot, wellRadius, params);

Assert(result.EndPressure, 1.7245e6, 'VerticalInjection_CO2_Pressure');
Assert(result.EndTemp, 156.08, 'VerticalInjection_CO2_Temp');
Assert(result.EndEnthalpy, 6.1802e5, 'VerticalInjection_CO2_Enthalpy');


% Horizontal well segment
T_res = T_e_initial + (params.dT_dz * abs(dz_total_v));
result = semi_analytic_well(result.EndPressure, result.EndTemp, T_res, dz_total_h, dr_total_h, m_dot, wellRadius, params);

Assert(result.EndPressure, 1.7238e6, 'HorizontalInjection_CO2_Pressure');
Assert(result.EndTemp, 212.746, 'HorizontalInjection_CO2_Temp');
Assert(result.EndEnthalpy, 6.755e5, 'HorizontalInjection_CO2_Enthalpy');

% High Flowrate
% Vertical well
m_dot = 150;
result = semi_analytic_well(P_f_initial, T_f_initial, T_e_initial, dz_total_v, dr_total_v, m_dot, wellRadius, params);

Assert(result.EndPressure, 5.1170e5, 'VerticalInjection_Highflow_CO2_Pressure');
Assert(result.EndTemp, 63.9786, 'VerticalInjection_Highflow_CO2_Temp');
Assert(result.EndEnthalpy, 5.3677e5, 'VerticalInjection_Highflow_CO2_Enthalpy');

% Small Radius
% Vertical well
m_dot = 0.1;
wellRadius = 0.02;
result = semi_analytic_well(P_f_initial, T_f_initial, T_e_initial, dz_total_v, dr_total_v, m_dot, wellRadius, params);

Assert(result.EndPressure, 1.1431e6, 'VerticalInjection_SmallRadius_CO2_Pressure');
Assert(result.EndTemp, 222.2246, 'VerticalInjection_SmallRadius_CO2_Temp');
Assert(result.EndEnthalpy, 6.8717e5, 'VerticalInjection_SmallRadius_CO2_Enthalpy');
