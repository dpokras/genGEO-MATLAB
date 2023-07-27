clear;
close all;
% Purpose of sheet is to generate a plot for the different system
% configurations (1 Lateral, 2 Laterals,...) for a variety of reservoir
% depths
%% Define System Parameters
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
num_Power_Plants = 2;
params.wellFieldType = 'Doublet'; 
params.time_years = 30;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';
% Select Fluid Type
fluid_All = ["Water","CO2"];
params.fluid = fluid_All(2);
%% Extract Depth and Corresponding Optimal Reservoir Lengths
C1_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction1_dTdz35_radius0.25.xlsx');
Res_Length_All = zeros(length(C1_Table(:,2)),4);
depth = C1_Table(:,1);
Res_Length_All(:,1) = C1_Table(:,2);
C2_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction2_dTdz35_radius0.25.xlsx');
Res_Length_All(:,2) = C2_Table(:,2);
C4_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz35_radius0.25.xlsx');
Res_Length_All(:,3) = C4_Table(:,2);
C8_Table = xlsread('data\Optimum_Res_Length_CO2_Conduction8_dTdz35_radius0.25.xlsx');
Res_Length_All(:,4) = C8_Table(:,2);
%% Select Variables in Plot
sys = ["Conduction1","Conduction2","Conduction4","Conduction8"];
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
    for j = 1:4
        params.system = sys(j);
        params.res_length = Res_Length_All(i,j);
        % Compute Power for each system type at some reservoir length
        result = total_analytic_system_optmdot(params);
        if sys(j) == "Conduction1"
            L = num_Power_Plants * (params.res_length * 1 + params.depth * 2);
        elseif sys(j) == "Conduction2"
            L = num_Power_Plants * (params.res_length * 2 + params.depth * 2);
        elseif sys(j) == "Conduction4"
            L = num_Power_Plants * (params.res_length * 4 + params.depth * 2);
        elseif sys(j) == "Conduction8"
            L = num_Power_Plants * (params.res_length * 8 + params.depth * 2);
        end
        W(i,j) = num_Power_Plants * result.W_net_IP / 1e6; % in MW
        Wm(i,j) = num_Power_Plants * result.W_net_IP / L; % in Watts/meter
        T_prod(i,j) = result.T_prod_surface_C;
        LCOE(i,j) = result.LCOE_greenfield * 1e6; %in 2019$/MWh
        m_dot(i,j) = result.m_dot_IP;
    end
end
%% Extract Data from Calculation and Input into .xlsx
Conduction_1_data = [depth,Res_Length_All(:,1),m_dot(:,1),W(:,1),Wm(:,1),T_prod(:,1),LCOE(:,1)/RF];
Conduction_2_data = [depth,Res_Length_All(:,2),m_dot(:,2),W(:,2),Wm(:,2),T_prod(:,2),LCOE(:,2)/RF];
Conduction_4_data = [depth,Res_Length_All(:,3),m_dot(:,3),W(:,3),Wm(:,3),T_prod(:,3),LCOE(:,3)/RF];
Conduction_8_data = [depth,Res_Length_All(:,4),m_dot(:,4),W(:,4),Wm(:,4),T_prod(:,4),LCOE(:,4)/RF];
toc;
T1 = table(Conduction_1_data);
writetable(T1,'data\Num_RWs_Conduction1.xlsx');
T2 = table(Conduction_2_data);
writetable(T2,'data\Num_RWs_Conduction2.xlsx');
T3 = table(Conduction_4_data);
writetable(T3,'data\Num_RWs_Conduction4.xlsx');
T8 = table(Conduction_8_data);
writetable(T8,'data\Num_RWs_Conduction8.xlsx');
