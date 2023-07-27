clear;
close all;
format short;
clc;
%%
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
params.system = 'Conduction4';
%params.wellCostType = 'Ideal';
params.wellCostType = 'Baseline';
params.dT_dz = 0.035;
params.k_rock = 2.1;
%d = 0.1016; %4 inch in meters
%d = 0.1143; %4.5 inch in meters
%d = 0.127; %5.00 inch in meters
%d = 0.13335; %5.25 inch in meters
%d = 0.1524; %6 inch in meters
%d = 0.2032; %8 inch in meters
%d = 0.244475; %9-5/8 inch in meters
d = 0.5; %base case
params.well_radius = d/2;

%% Define System Parameters
params.wellFieldType = 'Doublet'; 
params.time_years = 30;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';

%% Select Fluid Type
fluid_All = ["Water","CO2"];
%params.fluid = fluid_All(1);
params.fluid = fluid_All(2);
%% Initialize vectors
% depth = [7500];
depth = ((1:0.5:8)*1e3)';

res_length = zeros(1,length(depth));
LCOE = zeros(1,length(depth));

%% Find optimum reservoir length that minimizes LCOE at given depth
for i = 1:length(depth)
    params.depth = depth(i);
    params.res_length = 50;
    disp(strcat(['For depth of: ' num2str(params.depth) ]));
    result = total_analytic_system_optreslength(params);
    res_length(i) = result.res_length;
    LCOE(i) = result.LCOE_greenfield*1e6;
    disp(strcat(['// At a depth of: ' num2str(params.depth) ]));
    disp(strcat(['   Found a minimal LCOE of: ' num2str(result.LCOE_greenfield) ', at reservoir length: ' num2str(result.res_length) ]));
end
%% Output Excel File and Figure
depth = depth';
opt_res_length = res_length';
min_LCOE = LCOE';
T= table(depth, opt_res_length, min_LCOE);
if length(depth) > 1
    writetable(T,'data\k2_8.xlsx')
end
figure(1)
plot(res_length,depth','-k')
title('Optimum reservoir legths for specified depths using LCOE Ideal(CO_{2})')
xlabel('Resevoir Length [m]') % x-axis label
ylabel('Resevoir Depth [m]') % y-axis label
grid on