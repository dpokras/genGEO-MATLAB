clear;
close all;
clc;

% This sheet uses the optimal reservoir lengths determined for depths 1-8km
% for various therml gradients. Purpose of sheet is to find the corresponding
% power and LCOE and output it to a .xlsx file.
%% Define System Parameters
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
num_Power_Plants = 1;
params.wellFieldType = 'Doublet';
params.time_years = 30;
params.dT_dz = 0.035;
params.well_radius = 0.25;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';
params.system = 'Conduction4';
% Select Fluid Type
fluid_All = ["Water","CO2"];
params.fluid = fluid_All(2);
%% Extract Depth and Corresponding Optimal Reservoir Lengths
k2_1_Table = xlsread('data\Optimum_Res_Length_CO2_k2_1.xlsx');
Res_Length_All = zeros(length(k2_1_Table(:,2)),4);
depth = k2_1_Table(:,1);
Res_Length_All(:,1) = k2_1_Table(:,2);
k2_8_Table = xlsread('data\Optimum_Res_Length_CO2_k2_8.xlsx');
Res_Length_All(:,2) = k2_8_Table(:,2);
k3_5_Table = xlsread('data\Optimum_Res_Length_CO2_k3_5.xlsx');
Res_Length_All(:,3) = k3_5_Table(:,2);
%% Select Variables in Plot
k_T = [2.1,2.8,3.5];
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
    for j = 1:3 %3 is number of dTdz values
        params.k_rock = k_T(j);
        params.res_length = Res_Length_All(i,j);
        % Compute Outputs for each system type at some reservoir length
        result = total_analytic_system_optmdot(params);
        L = (params.res_length * 4 + params.depth * 2);
        W(i,j) =  result.W_net_IP / 1e6; % in MW
        Wm(i,j) = result.W_net_IP / L; % in Watts/meter
        T_prod(i,j) = result.T_prod_surface_C;
        LCOE(i,j) = result.LCOE_greenfield * 1e6; %in 2019$/MWh
        m_dot(i,j) = result.m_dot_IP;
    end
end
%% Extract Data from Calculation and Input into .xlsx
k_2_1_data = [depth,Res_Length_All(:,1),m_dot(:,1),W(:,1),Wm(:,1),T_prod(:,1),LCOE(:,1)];
k_2_8_data = [depth,Res_Length_All(:,2),m_dot(:,2),W(:,2),Wm(:,2),T_prod(:,2),LCOE(:,2)];
k_3_5_data = [depth,Res_Length_All(:,3),m_dot(:,3),W(:,3),Wm(:,3),T_prod(:,3),LCOE(:,3)];
T1 = table(k_2_1_data);
writetable(T1,'data\k_2_1.xlsx');
T2 = table(k_2_8_data);
writetable(T2,'data\k_2_8.xlsx');
T3 = table(k_3_5_data);
writetable(T3,'data\k_3_5.xlsx');
