% clear; clc; close all;

%%
times_yrs_init = 1:20:60;             % [0:.01:0.99 1:1:50];
times_yrs_full = 1:5:60;              % [0:.01:0.99 1:1:50];
times_months   = 0.1:0.1:1.1;         % [0:.01:0.99 1:1:50];

% parameters
fluid_id     = 2;
radius_id    = 2;
laterals     = 1;

working_fluids          = ["Water","CO2"];
well_radius_inch        = [12.25/2, 8.5/2];
well_radius_m           = well_radius_inch/39.37;

%% loading genGEO simulation parameters

params                  = SimulationParameters; 
params.optimizationMode = 'MaximizePower'; % MinimizeLCOE_Greenfield, MinimizeLCOE_Brownfield,  MaximizePower
params.system           = 'Conduction8';  % Porous, Conduction4
%params.thickness       = 100;
%params.transmissivity  = 1e-15 * 10^4.2;

%fixed inputs  ============================================================
params.depth            = 6000;
params.dT_dz            = 0.030;
params.T_surface_air_C  = 15;
params.T_surface_rock_C = 15;
params.dT_approach      = 7;

params.k_rock           = 3;            % thermal conductivity W/m-C 
params.C_rock           = 1000;         % Specific Heat J/kg-C
params.rho_rock         = 2650;         % Rock Density  kg/m3

% Variable inputs   =======================================================
params.fluid            = working_fluids(fluid_id);    
params.well_radius      = well_radius_m(radius_id);
params.res_length       = 1000;
params.m_dot_IP         = 15;

% Default inputs   ========================================================
params.silicaPrecip     = 'IgnoreSilica'; % PreventSilica, IgnoreSilica
%params.fieldType        = 'Brownfield'; % Greenfield, Brownfield
params.coolingMode      = 'Wet';      % Dry, Wet
params.orcFluid         = 'R245fa';   % R600a, R245fa
params.wellFieldType    = 'Doublet'; % Tungsten, Doublet, 5spot_SharedNeighbor
params.wellCostType     = 'Baseline'; %Ideal, Baseline

params.N_5spot          = 0;
params.F_OM             = 0.07;
params.discountRate     = 0.096;
params.Lifetime         = 25;
params.CapacityFactor   = 0.9;
params.costYear         = 2019;
%well_radius_vertical   = 0.25;
%well_radius_horizontal = 0.25;

%% Calculate the optimal mass flowrate and reservoir length (horizointal laterial length)
% for CO2-based AGS

optimal_res_len_co2 = zeros(0);

for j=1:length(init_res_len)
    for i=1:length(init_n_streams)
        
        params.time_years           = 30;
        params.res_length           = init_res_len(j);
        params.n_streams            = init_n_streams(i);
        params.fluid                = working_fluids(fluid_id);

        
        
        results_gengeo_yrs_optimal_co2 = total_analytic_system_optmdot(params);

        W_net_IP_MW     = results_gengeo_yrs_optimal_co2.W_net_IP/1e6;    
        Q_fluid_IP_MW   = results_gengeo_yrs_optimal_co2.Q_fluid_IP/1e6;
        %dP_surface_MPa  = results.dP_surface / 1e6;

        SpecificCapitalCost_brownfield_kW   = results_gengeo_yrs_optimal_co2.SpecificCapitalCost_brownfield * 1e3;
        LCOE_brownfield_MWh                 = results_gengeo_yrs_optimal_co2.LCOE_brownfield * 1e6;
        SpecificCapitalCost_greenfield_kW   = results_gengeo_yrs_optimal_co2.SpecificCapitalCost_greenfield * 1e3;
        LCOE_greenfield_MWh                 = results_gengeo_yrs_optimal_co2.LCOE_greenfield * 1e6;

        disp(strcat(['Power: ' num2str(W_net_IP_MW) ', LCOE: ' num2str(LCOE_brownfield_MWh)]));
        optimal_res_len_co2 = [optimal_res_len_co2; params.n_streams params.res_length   params.m_dot_IP W_net_IP_MW LCOE_brownfield_MWh results_gengeo_yrs_optimal_co2.m_dot_IP results_gengeo_yrs_optimal_co2.T_prod_surface_C];    
    % 
    end
end

%%
output_optimal_res_len_co2_n_streams            = optimal_res_len_co2(:,1);
output_optimal_res_len_co2_m_dot_ini           = optimal_res_len_co2(:,2);
output_optimal_res_len_co2_res_length          = optimal_res_len_co2(:,3);
output_optimal_res_len_co2_W_net_IP_MW         = optimal_res_len_co2(:,4);
output_optimal_res_len_co2_LCOE_brownfield_MWh = optimal_res_len_co2(:,5);
output_optimal_res_len_co2_m_dot_IP            = optimal_res_len_co2(:,6);
output_optimal_res_len_co2_T_prod_surface_C    = optimal_res_len_co2(:,7);

table_output_optimal_res_len_co2_time_yrs = table(output_optimal_res_len_co2_n_streams, ...
    output_optimal_res_len_co2_m_dot_ini, ...
    output_optimal_res_len_co2_res_length, ...
    output_optimal_res_len_co2_W_net_IP_MW,...
    output_optimal_res_len_co2_LCOE_brownfield_MWh,...
    output_optimal_res_len_co2_m_dot_IP,...
    output_optimal_res_len_co2_T_prod_surface_C);

writetable(table_output_optimal_res_len_co2_time_yrs,'C:\Users\pokra\Documents\ETH\Master Thesis\genGEO_matlab_AGS\yrs_output_optimal_res_len_co2_AGS_3.5kmdepth_4lateral_08inchdiam.xlsx')

% 
%% ========================= PLOT ResLength vs ProdTemp ===================
interval_streams = length(init_n_streams);
yr_1 = optimal_res_len_co2(1:interval_streams:end,:);
yr_2 = optimal_res_len_co2(2:interval_streams:end,:);
yr_3 = optimal_res_len_co2(3:interval_streams:end,:);


figure;
pos = get(gcf, 'Position');
set(gcf, 'Position',pos+[0 -500 0 500])

h = tiledlayout(3,1, 'TileSpacing', 'none', 'Padding', 'none');
nexttile

subplot(3,1,1);

plt_t01 = plot(yr_1(:,2),yr_1(:,5),'--ok','DisplayName','Prod Temp at Yr 01','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
hold on
plt_t21 = plot(yr_2(:,2),yr_2(:,5),':sk','DisplayName','Prod Temp at Yr 21', 'LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
plt_t41 = plot(yr_3(:,2),yr_3(:,5),'-k','Marker','diamond','DisplayName','Prod Temp at Yr 41', 'LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);

ylabel('LCOE [Eu/kWh]');
legend('Location','northwest');

grid on;
title(strcat({'Results for 3.5km deep conduction-based geothermal system using CO_2 with '...
    '4 horizontal lateral, 30 C/km gradient, and 08.5 Inch diameter at Year-1, 21,and  41.'}), 'FontSize', 10);


subplot(3,1,2);

plt_m01 = plot(yr_1(:,2),yr_1(:,6),'--ok','DisplayName','Mass flowrate at Yr 01','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
hold on
plt_m21 = plot(yr_2(:,2),yr_2(:,6),':sk','DisplayName','Mass flowrate  at Yr 21', 'LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
plt_m41 = plot(yr_3(:,2),yr_3(:,6),'-k','Marker','diamond','DisplayName','Mass flowrate  at Yr 41', 'LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);

ylabel('Mass Flowrate [kg/s]');
legend('Location','northwest');

grid on;

subplot(3,1,3);
%yyaxis left;
plt_p01 = plot(yr_1(:,2),yr_1(:,4),'--ok','DisplayName','Net Power at Yr 01','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
hold on
plt_p21 = plot(yr_2(:,2),yr_2(:,4),':sk','DisplayName','Net Power at Yr 21', 'LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
plt_p41 = plot(yr_3(:,2),yr_3(:,4),'-k','Marker','diamond','DisplayName','Net Power at Yr 41', 'LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);

ylabel('Power [MWe]');
xlabel( 'Horizontal Lateral Length [m]')

legend('Location','northwest');

grid on;

f = gcf;
exportgraphics(f,'C:\Users\pokra\Documents\ETH\Master Thesis\genGEO_matlab_AGS\output_optimal_res_len_co2_AGS_3.5kmdepth_4lateral_08inchdiameter.pdf','ContentType','vector')

%% Solving the power output for many years using the optimal value for the mass flow rate and lateral lenght
% in order to achieve the goal of 0.236 MWe after 40 years

opt_m_dot_solution = 35.62;
opt_res_length_solution = 5.5e3;

co2_gengeo_yrs_solution_output = zeros(0);

for i = 1:size(times_yrs_full, 2)
    params.time_years = times_yrs_full(i);
    params.m_dot_IP   = opt_m_dot_solution;
    params.res_length = opt_res_length_solution;
    params.fluid      = working_fluids(fluid_id);
    
    disp(strcat(['Trying year ' num2str(params.time_years)]));
    
    results_gengeo_yrs_solution_co2         = total_analytic_system_optmdot(params);
    
    W_net_IP_MW     = results_gengeo_yrs_solution_co2.W_net_IP / 1e6;    
    Q_fluid_IP_MW   = results_gengeo_yrs_solution_co2.Q_fluid_IP / 1e6;
    %dP_surface_MPa  = results.dP_surface / 1e6;
    
    SpecificCapitalCost_brownfield_kW   = results_gengeo_yrs_solution_co2.SpecificCapitalCost_brownfield * 1e3;
    LCOE_brownfield_MWh                 = results_gengeo_yrs_solution_co2.LCOE_brownfield * 1e6;
    SpecificCapitalCost_greenfield_kW   = results_gengeo_yrs_solution_co2.SpecificCapitalCost_greenfield * 1e3;
    LCOE_greenfield_MWh                 = results_gengeo_yrs_solution_co2.LCOE_greenfield * 1e6;
    
    disp(strcat(['Power: ' num2str(W_net_IP_MW) ', LCOE: ' num2str(LCOE_brownfield_MWh)]));
    co2_gengeo_yrs_solution_output = [co2_gengeo_yrs_solution_output; params.time_years W_net_IP_MW LCOE_brownfield_MWh results_gengeo_yrs_solution_co2.m_dot_IP results_gengeo_yrs_solution_co2.T_prod_surface_C];    
    
end


%% Solving the power output for first year (month-based) using the optimal value for the mass flow rate and lateral lenght

co2_gengeo_mnth_solution_output =zeros(0);

for i = 1:size(times_months, 2)
    params.time_years = times_months(i);
    params.m_dot_IP   = opt_m_dot_solution;
    params.res_length = opt_res_length_solution;
    params.fluid      = working_fluids(fluid_id);
    
    %disp(strcat(['Trying year ' num2str(params.time_years)]));
    
    results_gengeo_mnths_solution_co2         = total_analytic_system_optmdot(params);
    
    W_net_IP_MW     = results_gengeo_mnths_solution_co2.W_net_IP / 1e6;    
    Q_fluid_IP_MW   = results_gengeo_mnths_solution_co2.Q_fluid_IP / 1e6;
    %dP_surface_MPa  = results.dP_surface / 1e6;
    
    SpecificCapitalCost_brownfield_kW   = results_gengeo_mnths_solution_co2.SpecificCapitalCost_brownfield * 1e3;
    LCOE_brownfield_MWh                 = results_gengeo_mnths_solution_co2.LCOE_brownfield * 1e6;
    SpecificCapitalCost_greenfield_kW   = results_gengeo_mnths_solution_co2.SpecificCapitalCost_greenfield * 1e3;
    LCOE_greenfield_MWh                 = results_gengeo_mnths_solution_co2.LCOE_greenfield * 1e6;
    
    disp(strcat(['Power: ' num2str(W_net_IP_MW) ', LCOE: ' num2str(LCOE_brownfield_MWh)]));
    co2_gengeo_mnth_solution_output = [co2_gengeo_mnth_solution_output; params.time_years W_net_IP_MW LCOE_brownfield_MWh results_gengeo_mnths_solution_co2.m_dot_IP results_gengeo_mnths_solution_co2.T_prod_surface_C];    
    
end

%% Comparing the CO2 gengeo output to water on full year based

water_gengeo_yrs_solution_output=zeros(0);

for i = 1:size(times_yrs_full, 2)
    params.time_years   = times_yrs_full(i);
    params.fluid        = "Water";
    params.m_dot_IP     = opt_m_dot_solution;
    params.res_length   = opt_res_length_solution;
    
    disp(strcat(['Trying year ' num2str(params.time_years)]));
    
    results_gengeo_yrs_solution_water         = total_analytic_system_optmdot(params);
    
    W_net_IP_MW     = results_gengeo_yrs_solution_water.W_net_IP / 1e6;    
    Q_fluid_IP_MW   = results_gengeo_yrs_solution_water.Q_fluid_IP / 1e6;
    %dP_surface_MPa  = results.dP_surface / 1e6;
    
    SpecificCapitalCost_brownfield_kW   = results_gengeo_yrs_solution_water.SpecificCapitalCost_brownfield * 1e3;
    LCOE_brownfield_MWh                 = results_gengeo_yrs_solution_water.LCOE_brownfield * 1e6;
    SpecificCapitalCost_greenfield_kW   = results_gengeo_yrs_solution_water.SpecificCapitalCost_greenfield * 1e3;
    LCOE_greenfield_MWh                 = results_gengeo_yrs_solution_water.LCOE_greenfield * 1e6;
    
    disp(strcat(['Power: ' num2str(W_net_IP_MW) ', LCOE: ' num2str(LCOE_brownfield_MWh)]));
    water_gengeo_yrs_solution_output = [water_gengeo_yrs_solution_output; params.time_years W_net_IP_MW LCOE_brownfield_MWh results_gengeo_yrs_solution_water.m_dot_IP results_gengeo_yrs_solution_water.T_prod_surface_C];    
    
end

%% Comparing the CO2 gengeo output to water on first year (month  based)

water_gengeo_mnths_solution_output=zeros(0);

for i = 1:size(times_months, 2)
    params.time_years = times_months(i);
    params.fluid        = "Water";
    params.m_dot_IP     = opt_m_dot_solution;
    params.res_length   = opt_res_length_solution;
    
    disp(strcat(['Trying year ' num2str(params.time_years)]));
    
    results_gengeo_mnths_solution_water         = total_analytic_system_optmdot(params);
    
    W_net_IP_MW     = results_gengeo_mnths_solution_water.W_net_IP / 1e6;    
    Q_fluid_IP_MW   = results_gengeo_mnths_solution_water.Q_fluid_IP / 1e6;
    %dP_surface_MPa  = results.dP_surface / 1e6;
    
    SpecificCapitalCost_brownfield_kW   = results_gengeo_mnths_solution_water.SpecificCapitalCost_brownfield * 1e3;
    LCOE_brownfield_MWh                 = results_gengeo_mnths_solution_water.LCOE_brownfield * 1e6;
    SpecificCapitalCost_greenfield_kW   = results_gengeo_mnths_solution_water.SpecificCapitalCost_greenfield * 1e3;
    LCOE_greenfield_MWh                 = results_gengeo_mnths_solution_water.LCOE_greenfield * 1e6;
    
    disp(strcat(['Power: ' num2str(W_net_IP_MW) ', LCOE: ' num2str(LCOE_brownfield_MWh)]));
    water_gengeo_mnths_solution_output = [water_gengeo_mnths_solution_output; params.time_years W_net_IP_MW LCOE_brownfield_MWh results_gengeo_mnths_solution_water.m_dot_IP results_gengeo_mnths_solution_water.T_prod_surface_C];    
    
end

%% ========================PLOT===========================================

% plot the full-years-term
 
figure;
pos = get(gcf, 'Position');
set(gcf, 'Position',pos+[0 -500 0 500])

h = tiledlayout(3,1, 'TileSpacing', 'none', 'Padding', 'none');
nexttile

subplot(3,1,1);
plt_p01_co2 = plot(co2_gengeo_yrs_solution_output(:,1),co2_gengeo_yrs_solution_output(:,2),'-k','DisplayName','Net Power from CO_2','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
hold on
plt_p01_h2o = plot(water_gengeo_yrs_solution_output(:,1),water_gengeo_yrs_solution_output(:,2),'-r','DisplayName','Net Power from Water','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
 
ylabel('Power [MWe]');
legend('Location','best');

grid on;
title(strcat({'Results for 3.5km deep conduction-based geothermal system using CO_2 and '...
    ' Water with 4 horizontal lateral, 30 C/km gradient, and 08.5 Inch diameter over 60 Yrs.'}), 'FontSize', 10);

subplot(3,1,2);
plt_m01 = plot(co2_gengeo_yrs_solution_output(:,1),co2_gengeo_yrs_solution_output(:,4),'-k','DisplayName','Mass flowrate for CO_2','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
hold on
plt_m01 = plot(water_gengeo_yrs_solution_output(:,1),water_gengeo_yrs_solution_output(:,4),'-r','DisplayName','Mass flowrate for Water','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
grid on;

ylabel('Mass Flowrate [kg/s]');
legend('Location','best');

subplot(3,1,3);
plt_t01 = plot(co2_gengeo_yrs_solution_output(:,1),co2_gengeo_yrs_solution_output(:,5),'-k','DisplayName','Prod Temp for CO_2','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
hold on
plt_t01 = plot(water_gengeo_yrs_solution_output(:,1),water_gengeo_yrs_solution_output(:,5),'-r','DisplayName','Prod Temp for Water','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);

ylim([60,140])
ylabel('Production Temperature [C]');
legend('Location','best');

grid on

xlabel( 'Time [Years]')
%%
f = gcf;
exportgraphics(f,'C:\Users\pokra\Documents\ETH\Master Thesis\genGEO_matlab_AGS\power_plot_co2_water_many_yrs_3.5kmdepth_08inchdiameter.pdf','ContentType','vector')

%% plot the short-term (12 months)

figure;
pos = get(gcf, 'Position');
set(gcf, 'Position',pos+[0 -500 0 500])

h = tiledlayout(3,1, 'TileSpacing', 'none', 'Padding', 'none');
nexttile

subplot(3,1,1);
plt_p01_co2 = plot(co2_gengeo_mnth_solution_output(:,1),co2_gengeo_mnth_solution_output(:,2),'-k','DisplayName','Net Power from CO_2','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
hold on
plt_p01_h2o = plot(water_gengeo_mnths_solution_output(:,1),water_gengeo_mnths_solution_output(:,2),'-r','DisplayName','Net Power from Water','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);

ylabel('Power [MWe]');
legend('Location','best');

grid on;
title(strcat({'Results for 3.5km deep conduction-based geothermal system using CO_2 and '...
    ' Water with 4 horizontal lateral, 30 C/km gradient, and 08.5 Inch diameter over 1st Yr.'}), 'FontSize', 10);

subplot(3,1,2);
plt_m01 = plot(co2_gengeo_mnth_solution_output(:,1),co2_gengeo_mnth_solution_output(:,4),'-k','DisplayName','Mass flowrate for CO_2','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
hold on
plt_m01 = plot(water_gengeo_mnths_solution_output(:,1),water_gengeo_mnths_solution_output(:,4),'-r','DisplayName','Mass flowrate for Water','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
grid on;

ylabel('Mass Flowrate [kg/s]');
legend('Location','best');

subplot(3,1,3);
plt_t01 = plot(co2_gengeo_mnth_solution_output(:,1),co2_gengeo_mnth_solution_output(:,5),'-k','DisplayName','Prod Temp for CO_2','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);
hold on
plt_t01 = plot(water_gengeo_mnths_solution_output(:,1),water_gengeo_mnths_solution_output(:,5),'-r','DisplayName','Prod Temp for Water','LineWidth',2, 'MarkerFaceColor',[0 0 0], 'MarkerSize',5);

ylim([60,140])
ylabel('Production Temperature [C]');
legend('Location','best');

grid on
xlabel( 'Time [Years]')

f = gcf;
exportgraphics(f,'gold\final\power_plot_co2_water_many_months_3.5kmdepth_08inchdiameter.pdf','ContentType','vector')


