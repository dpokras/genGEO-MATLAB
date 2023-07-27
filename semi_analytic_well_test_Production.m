clear all;
close all;

params = SimulationParameters;
params.dT_dz = 0.035;
params.time_years = 10;
params.fluid = 'Water';

P_f_initial = 25e6;
T_f_initial = 97;
T_e_initial = 102.5;
dz_total = 2500;
dr_total = 0;
m_dot = 136;
wellRadius = 0.205;

% Vertical well segment
result = semi_analytic_well(P_f_initial, T_f_initial, T_e_initial, dz_total, dr_total, m_dot, wellRadius, params);

Assert(result.EndPressure, 1.2372e6, 'Production_Water_Pressure');
Assert(result.EndTemp, 95.0133, 'Production_Water_Temp');
Assert(result.EndEnthalpy, 3.9902e5, 'Production_Water_Enthalpy');

% Now CO2
params.fluid = 'CO2';
result = semi_analytic_well(P_f_initial, T_f_initial, T_e_initial, dz_total, dr_total, m_dot, wellRadius, params);

Assert(result.EndPressure, 1.1674e7, 'Production_CO2_Pressure');
Assert(result.EndTemp, 55.9783, 'Production_CO2_Temp');
Assert(result.EndEnthalpy, 3.7225e5, 'Production_CO2_Enthalpy');