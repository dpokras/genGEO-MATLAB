clear;
close all;
clc;

% This sheet uses the optimal reservoir lengths determined for depths 1-8km
% for various therml gradients. Purpose of sheet is to find the corresponding
% power and LCOE and output it to a .xlsx file.
%% Define System Parameters
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
num_Power_Plants = 2;
params.wellFieldType = 'Doublet'; 
params.time_years = 30;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';
params.system = 'Conduction4';
% Select Fluid Type
fluid_All = ["Water","CO2"];
params.fluid = fluid_All(2);

%% Extract Depth and Corresponding Optimal Reservoir Lengths
C1_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz25_radius0.25.xlsx');
Res_Length_All = zeros(length(C1_Table(:,2)),4);
depth = C1_Table(:,1);
Res_Length_All(:,1) = C1_Table(:,2);
C2_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz30_radius0.25.xlsx');
Res_Length_All(:,2) = C2_Table(:,2);
C4_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz35_radius0.25.xlsx');
Res_Length_All(:,3) = C4_Table(:,2);
C8_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz40_radius0.25.xlsx');
Res_Length_All(:,4) = C8_Table(:,2);

%% Extract Depth and Corresponding Optimal Reservoir Lengths
% C1_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz25_radius0.25.xlsx');
% Res_Length_All = zeros(length(C1_Table(:,2)),4);
% depth = C1_Table(:,1);
% Res_Length_All(:,1) = C1_Table(:,2);
% C2_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz30_radius0.25.xlsx');
% Res_Length_All(:,2) = C2_Table(:,2);
% C4_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz35_radius0.25.xlsx');
% Res_Length_All(:,3) = C4_Table(:,2);
% C8_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz40_radius0.25.xlsx');
% Res_Length_All(:,4) = C8_Table(:,2);
%% Select Variables in Plot
gradients = [0.025,0.030,0.035,0.040];
%% Initialize Vectors
W = zeros(length(depth),4);
Wm = W;
T_prod = W;
LCOE = W;
m_dot = W;

%% Calculate LCOE, Power, Specific Power, Production Temperature of System for Each Depth
tic;
for i = 1:length(depth) %for each resevoir length
    params.depth = depth(i);
    for j = 1:4 %4 is number of dTdz values
        params.dT_dz = gradients(j);
        params.res_length = Res_Length_All(i,j);
        % Compute Outputs for each system type at some reservoir length
        result = total_analytic_system_optmdot(params);
        L = num_Power_Plants * (params.res_length * 4 + params.depth * 2);
        W(i,j) = num_Power_Plants * result.W_net_IP / 1e6; % in MW
        Wm(i,j) = num_Power_Plants * result.W_net_IP / L; % in Watts/meter
        T_prod(i,j) = result.T_prod_surface_C;
        LCOE(i,j) = result.LCOE_greenfield * 1e6; %in 2019$/MWh
        m_dot(i,j) = result.m_dot_IP;
    end
end
%% Extract Data from Calculation and Input into .xlsx
dTdz_25_data = [depth,Res_Length_All(:,1),m_dot(:,1),W(:,1),Wm(:,1),T_prod(:,1),LCOE(:,1)];
dTdz_30_data = [depth,Res_Length_All(:,2),m_dot(:,2),W(:,2),Wm(:,2),T_prod(:,2),LCOE(:,2)];
dTdz_35_data = [depth,Res_Length_All(:,3),m_dot(:,3),W(:,3),Wm(:,3),T_prod(:,3),LCOE(:,3)];
dTdz_40_data = [depth,Res_Length_All(:,4),m_dot(:,4),W(:,4),Wm(:,4),T_prod(:,4),LCOE(:,4)];
toc;
T1 = table(dTdz_25_data);
writetable(T1,'data\dTdz_25.xlsx');
T2 = table(dTdz_30_data);
writetable(T2,'data\dTdz_30.xlsx');
T3 = table(dTdz_35_data);
writetable(T3,'data\dTdz_35.xlsx');
T8 = table(dTdz_40_data);
writetable(T8,'data\dTdz_40.xlsx');