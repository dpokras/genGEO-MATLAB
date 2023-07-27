clear; clc; close all;

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
params.optimizationMode = 'MinimizeLCOE_Greenfield'; %'MaximizePower', MinimizeLCOE_Brownfield,  MaximizePower
params.system           = 'Conduction8';  % Porous, Conduction4
%params.thickness       = 100;
%params.transmissivity  = 1e-15 * 10^4.2;

%fixed inputs  ============================================================
% params.depth            = 6000;
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
% params.res_length       = 1000;
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

depths = 2000:500:8000;
res_lens  = 500:500:6000;
n_streams = 1:5:21;


output = zeros(0);

for i=1:length(depths)
    for j=1:length(res_lens)
        for k=1:length(n_streams)
            
            params.time_years           = 30;
            params.depth                = depths(i)
            params.res_length           = res_lens(j);
            params.n_streams            = n_streams(k)
            params.fluid                = working_fluids(fluid_id);
            
            disp(strcat(['depth: ' num2str(depths(i)) ' res_length: ' num2str(res_lens(j)) ' no streams: ' num2str(n_streams(k))]))
            results = total_analytic_system_optmdot(params);
    
    
    
            W_net_IP_MW     = results.W_net_IP/1e6;    
            Q_fluid_IP_MW   = results.Q_fluid_IP/1e6;
            %dP_surface_MPa  = results.dP_surface / 1e6;
    
            SpecificCapitalCost_brownfield_kW   = results.SpecificCapitalCost_brownfield * 1e3;
            LCOE_brownfield_MWh                 = results.LCOE_brownfield * 1e6;
            SpecificCapitalCost_greenfield_kW   = results.SpecificCapitalCost_greenfield * 1e3;
            LCOE_greenfield_MWh                 = results.LCOE_greenfield * 1e6;
    
            disp(strcat(['Power: ' num2str(W_net_IP_MW) ', LCOE: ' num2str(LCOE_brownfield_MWh)]));
            output = [output; params.n_streams params.depth params.res_length W_net_IP_MW LCOE_brownfield_MWh results.m_dot_IP results.T_prod_surface_C];
        end
    end
end

%%
output_n_streams            = output(:,1);
output_res_length            = output(:,2);
output_W_net_IP_MW          = output(:,4);
output_LCOE_brownfield_MWh  = output(:,5);
output_m_dot_IP             = output(:,6);
output_T_prod_surface_C     = output(:,7);

output_table = table( ...
    output_n_streams, ...
    output_res_length, ...
    output_W_net_IP_MW,...
    output_LCOE_brownfield_MWh,...
    output_m_dot_IP,...
    output_T_prod_surface_C);

writetable(output_table,'C:\Users\pokra\Documents\ETH\Master Thesis\genGEO_matlab_AGS\yrs_output_optimal_res_len_co2_AGS_3.5kmdepth_4lateral_08inchdiam.xlsx')

% 
%% ========================= PLOT ResLength vs ProdTemp ===================
interval_streams = length(init_n_streams);
X = output(1:interval_streams:end,:);
Y = output(2:interval_streams:end,:);
Z = output(3:interval_streams:end,:);


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
