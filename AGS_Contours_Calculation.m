clear;
close all;
clc;
tic;

% First find the optimal m-dot. Then, use m to find Tprod
%% Define parameters for base case
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
params.system = 'Conduction4';

%% Define System Parameters
params.wellFieldType = 'Doublet';
params.time_years = 30;
params.wellCostType = 'Ideal';
%params.wellCostType = 'Baseline';
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';
fluid_All = ["Water","CO2"];
% Select Fluid Type
%params.fluid = fluid_All(1);
params.fluid = fluid_All(2);
%d = 0.1016; %4 inch in meters
%d = 0.1524; %6 inch in meters
%d = 0.2032; %8 inch in meters
%d = 0.244475; %9-5/8 inch in meters
d = 0.5; %base case
params.well_radius = d/2;
%% Select Extent of Depth and Reservoir Length
% Total Reservoir Length
res_length = linspace(0,30000,61);
% Total Reservoir Depth
total_depth = 8000;
depth = linspace(1000,total_depth,15);
%% Initialize Vectors
W_net_IP_kW = zeros(1, length(depth) * length(res_length));
LCOE_MWh = zeros(1, length(depth) * length(res_length));
W_specific = zeros(1, length(depth) * length(res_length));
T_prod = zeros(1, length(depth) * length(res_length));
counter = 0; % counter for indexing values in W_net_IP
res_vector = zeros(1, length(depth) * length(res_length));
depth_vector = zeros(1, length(depth) * length(res_length));
%% Find Power, Specific Power, LCOE, Production Temperature for all depths and resevoir lengths
for i = 1:length(depth)
    for j = 1:length(res_length)
        counter = counter + 1;
        depth_vector(counter) = depth(i);
        res_vector(counter) = res_length(j);
        if params.system == "Conduction4"
            L = 2*depth(i) + 4*res_length(j); % NOTE: DO not change unless something other than Conduction4 is used
        else
            throw(MException('AGS_Contours_Calculation:NotImplemented','AGS_Contours_Calculation:Redefine params.system'));
        end
            %% Find m-dot at the res length and depth
        params.depth = depth(i);
        params.res_length = res_length(j);
        disp(strcat(['For depth of: ' num2str(params.depth) ' and a reservoir length of: ' num2str(params.res_length) ]));
        result = total_analytic_system_optmdot(params);
        %% Use Optimal flow rate to find Production Temperature, Power Generated, Specific Power, LCOE
        params.m_dot_IP = result.m_dot_IP;
        W_net_IP_kW(counter) = result.W_net_IP / 1e3; %system of 1 with kWe
        LCOE_MWh(counter) = result.LCOE_greenfield * 1e6; % Do not change regardless of system
        W_specific(counter) = result.W_net_IP / L; %in We/m
        T_prod(counter) = result.T_prod_surface_C;
        disp(strcat(['// At a depth of: ' num2str(params.depth) ' and reservoir length of: ' num2str(params.res_length) ]));
        disp(strcat(['   Found a minimal LCOE of: ' num2str(LCOE_MWh(counter)) ', Power of: ' num2str(W_net_IP_kW(counter)) ]));
        disp(strcat(['   a specific power of: ' num2str(W_specific(counter)) ', Prod Temp. of: ' num2str(T_prod(counter)) ]));
    end
end
%% Write Tables
% For LCOE
T1 = table(depth_vector', res_vector', LCOE_MWh');
writetable(T1,'data\Contour_CO2_LCOE_Conduction4_dTdz35_r25_Ideal-Feb2021.xlsx')

% For Power
T2 = table(depth_vector', res_vector', W_net_IP_kW');
writetable(T2,'data\Contour_CO2_Power_Conduction4_dTdz35_r25_Ideal-Feb2021.xlsx')

% For Specific Power
T3 = table(depth_vector', res_vector', W_specific');
writetable(T3,'data\Contour_CO2_Specific_Power_Conduction4_dTdz35_r25_Ideal-Feb2021.xlsx')

% For Production Temperature
T4 = table(depth_vector', res_vector', T_prod');
writetable(T4,'data\Contour_CO2_Production_Temp_Conduction4_dTdz35_r25_Ideal-Feb2021.xlsx')

toc;
