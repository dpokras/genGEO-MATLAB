%% Financial Aspects
F_OM = 0.045;
d = 0.096;
n = 25;
CF = 0.85;
CRF = (d*(1+d)^n)/((1+d)^n-1);
RF = (CRF + F_OM) / (CF * 8760) * 1e6
depth = linspace(1,8,15)';
res = linspace(0,30,61);

LCOE_Baseline = xlsread('data\Contour_CO2_LCOE_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
LCOE_Ideal = xlsread('data\Contour_CO2_LCOE_Conduction4_dTdz35_r25_Ideal-Feb2021.xlsx');
Power = xlsread('data\Contour_CO2_Power_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
Power = Power(:,3)/2;
Specific_Power = xlsread('data\Contour_CO2_Specific_Power_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
Specific_Power = Specific_Power(:,3);
Temp = xlsread('data\Contour_CO2_Production_Temp_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
Temp = Temp(:,3);
depth_baseline = LCOE_Baseline(:,1)/1000;
depth_ideal = LCOE_Ideal(:,1)/1000;
res_baseline = LCOE_Baseline(:,2)/1000;
res_ideal = LCOE_Ideal(:,2)/1000;

SpCC_Baseline = LCOE_Baseline(:,3)/RF;
SpCC_Ideal = LCOE_Ideal(:,3)/RF;

value_table = zeros(length(depth)+1,length(res)+1);
value_table(2:end,1) = depth;
value_table(1,2:end) = res;

r = 1;
r_end = 61;
for i = 1:15
    z_values = Temp(r:r_end)';
    value_table(i+1,2:end) = z_values;
    r = r_end+1;
    r_end = r_end + 61;
end
value_table = array2table(value_table)
writetable(value_table,'data\Figure4_Data_Temp.xlsx')

% Power Data
CO2_Power = xlsread('data\Contour_CO2_Power_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
CO2_Power = CO2_Power(:,3)/2;
% Specific Power Data
CO2_Specific_Power = xlsread('data\Contour_CO2_Specific_Power_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
CO2_Specific_Power = CO2_Specific_Power(:,3);
% Production Temperature Data
CO2_Temp = xlsread('data\Contour_CO2_Production_Temp_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
CO2_Temp = CO2_Temp(:,3);
% Optimum Reservoir Length Data
CO2_Res_Length_Baseline = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz35_radius0.25.xlsx');
CO2_Res_Length_Ideal = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz35_radius25_ideal.xlsx');
opt_res_length_co2_baseline = CO2_Res_Length_Baseline(:,2)/1000;
opt_res_length_co2_ideal = CO2_Res_Length_Ideal(:,2)/1000;
depth_co2_baseline = CO2_Res_Length_Baseline(:,1)/1000;
depth_co2_ideal = CO2_Res_Length_Ideal(:,1)/1000;