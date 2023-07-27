clear;
close all;
clc;
tic;

% First find the optimal m-dot. Then, use m to find Tprod
%% Define parameters for base case
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
%% LW
params.res_length = 5000; % base case

%% DIAMETER
d = 0.244475; %9-5/8 inch in meters
%d = 0.5; %base case
params.well_radius = d/2;
%% DEPTH
params.depth = 8000;
%% ALL OTHER PARAMETERS
params.system = 'Conduction1';
params.wellFieldType = 'Doublet'; 
params.time_years = 30;
params.dT_dz = 0.035;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';
fluid_All = ["Water","CO2"];
% Select Fluid Type
params.fluid = fluid_All(2);
%% Calculate Outputs for Baseline
params.wellCostType = 'Baseline';
result = total_analytic_system_optmdot(params);
% output_baseline = [params.depth, params.res_length, result.m_dot_IP, result.W_net_IP/1e6, result.T_prod_surface_C, result.W_net_IP/(2*params.depth + 4*params.res_length), result.CapitalCost.C_greenfield,result.SpecificCapitalCost_greenfield];
% spcc = output_baseline(end)
% power = result.W_net_IP/1e6
% 
% Cost Output
cc_prod_baseline = result.CapitalCost.C_wells_production/1e6;
cc_inj_baseline = result.CapitalCost.C_wells_injection/1e6;
cc_drilling_baseline = cc_inj_baseline + cc_prod_baseline
cc_wellfield_baseline = result.CapitalCost.C_wellfield/1e6
cc_stimulation_baseline = result.CapitalCost.C_stimulation/1e6
cc_exploration_baseline = result.CapitalCost.C_exploration/1e6
cc_gatheringsystem_baseline = result.CapitalCost.C_gatheringSystem/1e6
cc_surfacePlant_baseline = result.CapitalCost.C_surfacePlant/1e6
cc_total_baseline = result.CapitalCost.C_greenfield/1e6
cc_total_verify_baseline = cc_drilling_baseline + cc_wellfield_baseline + cc_stimulation_baseline + cc_exploration_baseline + cc_gatheringsystem_baseline + cc_surfacePlant_baseline
cost_baseline = [cc_drilling_baseline, cc_wellfield_baseline,cc_stimulation_baseline, cc_exploration_baseline, cc_gatheringsystem_baseline, cc_surfacePlant_baseline, cc_total_baseline];
% 
% 
%% Calculate Outputs for Ideal
params.wellCostType = 'Ideal';
result = total_analytic_system_optmdot(params);
% output_ideal = [params.depth, params.res_length, result.m_dot_IP, result.W_net_IP/1e6, result.T_prod_surface_C, result.W_net_IP/(2*params.depth + 4*params.res_length), result.CapitalCost.C_greenfield, result.SpecificCapitalCost_greenfield]
% 
% %% Cost Output
cc_prod_ideal = result.CapitalCost.C_wells_production/1e6;
cc_inj_ideal = result.CapitalCost.C_wells_injection/1e6;
cc_drilling_ideal = cc_inj_ideal + cc_prod_ideal
cc_wellfield_ideal = result.CapitalCost.C_wellfield/1e6
cc_stimulation_ideal = result.CapitalCost.C_stimulation/1e6
cc_exploration_ideal = result.CapitalCost.C_exploration/1e6
cc_gatheringsystem_ideal = result.CapitalCost.C_gatheringSystem/1e6
cc_surfacePlant_ideal = result.CapitalCost.C_surfacePlant/1e6
cc_total_ideal = result.CapitalCost.C_greenfield/1e6
cc_total_verify_ideal = cc_drilling_ideal + cc_wellfield_ideal + cc_stimulation_ideal + cc_exploration_ideal + cc_gatheringsystem_ideal + cc_surfacePlant_ideal

cost_ideal = [cc_drilling_ideal, cc_wellfield_ideal,cc_stimulation_ideal, cc_exploration_ideal, cc_gatheringsystem_ideal, cc_surfacePlant_ideal, cc_total_ideal];

cost_table = [cost_baseline, 0, cost_ideal];
cost_table = array2table(cost_table)
writetable(cost_table,'data\Figure10_Pie_Chart.xlsx')
% percent_drilling_baseline = (cc_drilling_baseline)*100 / cc_total_baseline
% percent_drilling_ideal = (cc_drilling_ideal)*100 / cc_total_ideal

%% Pie Charts
figure(1)
subplot(2,1,1);
cc_drilling_baseline = 103.2;
cc_wellfield = 1.5;
cc_gatheringsystem = 1;
cc_surfacePlant = 3.1;
labels = {'Well Drilling & Completion','Wellfield','Gathering System','Surface Power Plant'};
baseline = [cc_drilling_baseline, cc_wellfield, cc_gatheringsystem, cc_surfacePlant];
label_pie = {'103.2M$ (95%)','1.5M$ (1%)','1.0M$ (1%)','3.1M$ (3%)'};

explode = [1,1,1,1];
pie(baseline,label_pie)
get(gca, 'position');
ax1 = gca;
ax1.View = [10 90];
set(ax1, 'Position', [0.13 0.5 0.8 0.45],'fontsize',11);
legend(labels)
title('Baseline Capital Costs')

subplot(2,1,2);
cc_drilling_ideal = 28;
labels = {'Well Drilling & Completion','Wellfield','Gathering System','Surface Power Plant'};
ideal = [cc_drilling_ideal, cc_wellfield, cc_gatheringsystem, cc_surfacePlant];
label_pie = {'28.0M$ (84%)','1.5M$ (4%)','1.0M$ (3%)','3.1M$ (9%)'};
pie(ideal, label_pie)
ax2 = gca;
get(ax2, 'position');
ax2.View = [10 90];
set(ax2, 'Position', [0.13 0.01 0.8 0.45],'fontsize',11);
legend(labels)
title('Ideal Capital Costs')

caption = {'{\bfProperties}'...
    ,'8.0 km Vertical Well Depth'...
    ,'5.0 km Lateral Well Length'...
    ,'Massflow to Minimize Cost'...
    ,'CO_{2} Working Fluid'...
    ,'35 �C/km Temp Gradient'...
    ,'15 �C Surface Temp'...
    ,'9 5/8-inch Well Diameter'...
    ,'One Horizontal Laterals'...
    ,strcat(['Year-' num2str(30) ' Depletion Values'])...
    ,'Malek|Adams|Rossi|Schiegg|Saar (2021)'};
annot = annotation('textbox','String',caption);
annot.Position = [0.65 0.69 0.4 0.3];
annot.FitBoxToText = 'on';
annot.BackgroundColor = 'White';
annot.FontSize = 11;

saveas(figure(1),'images\PieCharts.fig');


