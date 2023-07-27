clear;
close all;
clc;
tic;

%% Define parameters for base case
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
num_Power_Plants = 1;
params.res_length = 5000; % in meters
params.depth = 3500; % in meters
params.system = 'Conduction4'; % Horizontal pipe section consists of 4 laterals
params.wellFieldType = 'Doublet';    
L = (4*params.res_length + 2*params.depth)*num_Power_Plants;

%% Define System Parameters
m_dot_IP_max = 100; % Input maximum flowrate
params.time_years = 30;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';
fluid_All = ["Water","CO2"];

%% Initialize Vectors
W_net_IP_MW = zeros(2,m_dot_IP_max+1); % Initialize power generated vector for water and CO2
specific_w = zeros(2,m_dot_IP_max+1); % Initialize specific power vector for water and CO2
% P_difference = zeros(3,m_dot_IP_max+1); % Initialize dP vector
%                                         % Row 1 for dP_pump water, 
%                                         % Row 2 for dP_turbine Co2
%                                         % Row 3 for dP_pump Co2
T_prod = zeros(2,m_dot_IP_max+1); % Initialize production temperature vector for water and CO2
dp_co2 = zeros(1,m_dot_IP_max+1);
SpCC_greenfield = zeros(2,m_dot_IP_max+1); %Initialize SpCC Greenfield vecor for water and CO2
m = linspace(0,m_dot_IP_max,m_dot_IP_max+1);
m = fliplr(m);
m_dot_opt = zeros(2,1);
W_max = zeros(2,1);
SpCC_min = zeros(2,1);
temp_prod = zeros(2,1);
specific_w_max = zeros(2,1);
%% For Each Fluid Type and for many m-dots, find Power, Specific Power, SpCC, Relevant Pressures 

for j = 1:2 % for each fluid type
    params.fluid = fluid_All(j);
    for i = 1:m_dot_IP_max+1
            params.m_dot_IP = m(i);
            if params.m_dot_IP == 0
                result.dP_turbine = 0;
                result.W_net_IP = 0;
                result.T_prod_surface_C = params.T_surface_rock_C;
            else    
                if params.fluid == "Water"
                    result = total_analytic_system_water(params);
                elseif params.fluid == "CO2"
                % CODE below uses total analytic system CO2 for some m-dot to find values
                    try
                     result = total_analytic_system_co2(params);
                    catch MException('CapitalCost_SurfacePlant:PositiveCondenserHeat','Positive Condenser Heat')
                        result.W_net_IP = 0;
                        m_interp = m(1:i-1);
                        T_interp = T_prod(j,1:i-1);
                        result.T_prod_surface_C = interp1(m_interp,T_interp,m(i),'linear','extrap');
                        P_interp =  P_difference(j,1:i-1);
                        result.dP_turbine = interp1(m_interp,P_interp,m(i),'linear','extrap');
                    end
                end
            end
            if params.fluid == "Water"
                P_difference(1,i) = result.dP_pump / 1e6; % Use dP_Pump if fluid is water
                W_net_IP_MW(1,i) = num_Power_Plants * result.W_net_IP / 1e6;
                T_prod(1,i) = result.T_prod_surface_C;
                SpCC_greenfield(1,i) = result.SpecificCapitalCost_greenfield;
                specific_w(1,i) = W_net_IP_MW(1,i)*1e6 / L;
            end
            if params.fluid == "CO2"
                P_difference(2,i) = result.dP_turbine / 1e6;
                P_difference(3,i) = result.dP_pump / 1e6;
                W_net_IP_MW(2,i) = num_Power_Plants * result.W_net_IP / 1e6;
                T_prod(2,i) = result.T_prod_surface_C;
                SpCC_greenfield(2,i) = result.SpecificCapitalCost_greenfield;
                specific_w(2,i) = W_net_IP_MW(2,i)*1e6 / L;
            end
    end
    
%% Determine Maximum Power Output and corresponding m-dot
    params.fluid = fluid_All(j);
    % Determine Optimum m-dot
    result = total_analytic_system_optmdot(params);
    m_dot_opt(j) = result.m_dot_IP;
    SpCC_min(j) = result.SpecificCapitalCost_greenfield;
    temp_prod(j) = result.T_prod_surface_C;
    specific_w_max(j) = result.W_net_IP / L;
    W_max(j) = num_Power_Plants * result.W_net_IP / 1e6;

%% Plot for Power Generated and Production Temperature
    figure(1)
    subplot(2,2,2)
    xlim([0,m_dot_IP_max]) % Specifying limits on x-axis
    xlabel('Mass Flowrate $\dot{m}$ [kg s$^{-1}$]', 'FontSize',11,'interpreter','latex') % x-axis label
    yyaxis left
    water = plot(NaN,NaN,'-k','linewidth',1.4);
    hold on
    co2 = plot(NaN,NaN,'--k','linewidth',1.4);
    optimal = plot(NaN,NaN,'pentagram','MarkerFaceColor','#000000','MarkerEdgeColor','#000000', 'MarkerSize',5);
    if params.fluid == "Water"
        T_water = plot(m,T_prod(1,:),'-','Color','#000000','linewidth',1.4);
        hold on 
        T_water_max = plot(m_dot_opt(1), temp_prod(1),'pentagram','MarkerFaceColor','#000000', 'MarkerEdgeColor','#000000','MarkerSize',5);
    end
    ylabel('Production Temperature [\circ{}C]') % y-axis label
    ylim([0,120]);
    hold on
    if params.fluid == "CO2"
        T_co2 = plot(m,T_prod(2,:),'--','Color','#000000','linewidth',1.4);
        hold on 
        T_co2_max = plot(m_dot_opt(2), temp_prod(2),'pentagram','MarkerFaceColor','#000000', 'MarkerEdgeColor','#000000','MarkerSize',5);
    end
    yyaxis right
     if params.fluid == "Water"
        W_water = plot(m,W_net_IP_MW(1,:),'-','Color','#4e5fa2','linewidth',1.4);
        opt_water = plot(m_dot_opt(1),W_max(1),'pentagram','MarkerFaceColor','#4e5fa2','MarkerEdgeColor','#4e5fa2', 'MarkerSize',5);
    end
    ylabel('Electric Power [MW_{e}]', 'FontSize',11); % y-axis label
    grid on
    hold on
    if params.fluid == "CO2"
        W_co2 = plot(m,W_net_IP_MW(2,:),'--', 'Color','#4e5fa2','linewidth',1.4);
        opt_co2 = plot(m_dot_opt(2),W_max(2),'pentagram','MarkerFaceColor','#4e5fa2', 'MarkerEdgeColor','#4e5fa2','MarkerSize',5);
    end
    ylim([0,0.4]) % Specifying limits on y-axis
    hold on
    if params.fluid == "CO2"
        lgd = legend([water co2 optimal], 'Water','CO_{2}', 'Optimal Mass Flowrate');
        lgd.Title.String = 'Legend';
        lgd.Title.FontSize = 11;
        lgd.FontSize = 11;
        lgd.Position = [0.15 0.6 0.3 0.3]; %[110 500 150 90];
    end
    axt = gca;
    axt.YAxis(1).Color = '#000000';
    axt.YAxis(2).Color = '#4e5fa2';
    set(axt,'FontSize',11, 'gridcolor','k');
    grid on
%% Plot For SpCC and Specific Power
    subplot(2,2,4)
    SpCC_greenfield(1,end) = SpCC_greenfield(1,1)^100; % To make SpCC @ m = 0 very high
    SpCC_greenfield(2,end)= SpCC_greenfield(2,1)^100; % To make SpCC @ m = 0 very high
    yyaxis left
    if params.fluid == "Water"
        SpCC_w = plot(m,SpCC_greenfield(1,:),'-k','linewidth',1.4);
        hold on
        SpCC_w_min = plot(m_dot_opt(1),SpCC_min(1),'pentagram','MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',5);
    end
    xlim([0,m_dot_IP_max]) % Specifying limits on x-axis
    ylim([0, 1500]) % Specifying limits on y-axis so that it doesn't go to infinity
    xlabel('Mass Flowrate $\dot{m}$ [kg s$^{-1}$]', 'FontSize',11,'interpreter','latex') % x-axis label
    ylabel('SpCC [$ W_{e}^{-1}]', 'FontSize',11) % y-axis label
    hold on
    if params.fluid == "CO2"
        SpCC_c = plot(m,SpCC_greenfield(2,:),'--k','linewidth',1.4);
        hold on
        SpCC_c_min = plot(m_dot_opt(2),SpCC_min(2),'pentagram','MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',5);
    end
    yyaxis right
    if params.fluid == "Water"
        spec_W_w = plot(m,specific_w(1,:),'-','Color','#4e5fa2','linewidth',1.4);
        hold on 
        spec_W_w_max = plot(m_dot_opt(1), specific_w_max(1),'pentagram','MarkerFaceColor','#4e5fa2', 'MarkerEdgeColor','#4e5fa2','MarkerSize',5);
    end
    if params.fluid == "CO2"
        spec_W_c = plot(m,specific_w(2,:),'--','Color','#4e5fa2','linewidth',1.4);
        hold on
        spec_W_c_max = plot(m_dot_opt(2), specific_w_max(2),'pentagram','MarkerFaceColor','#4e5fa2', 'MarkerEdgeColor','#4e5fa2','MarkerSize',5);
    end
    ylabel('Specific Electric Power [W_{e} m^{-1}]', 'FontSize',11) % y-axis label
    ylim([0,15]);
    hold on
    set(gca,'FontSize',11, 'gridcolor','k');
    axb = gca;
    axb.YAxis(1).Color = 'k';
    axb.YAxis(2).Color = '#4e5fa2';   
    grid on
end
set(axt, 'Position', [0.37 0.6 0.55 0.37])
set(axb, 'Position', [0.37 0.10 0.55 0.37])
caption = {'{\bfProperties}'...
    ,'3.5 km Vertical Well Depth'...
    ,'5.0 km Lateral Well Length'...
    ,'35 �C/km Temp Gradient'...
    ,'15 �C Surface Temp'...
    ,'0.5 m Well Diameter'...
    ,'Four Horizontal Laterals'...
    ,strcat(['Year-' num2str(30) ' Depletion Values'])...
    ,'Malek|Adams|Rossi|Schiegg|Saar (2021)'};
annot = annotation('textbox','String',caption);
annot.Position = [0.03 0.67 0.3 0.3];
annot.FitBoxToText = 'on';
annot.BackgroundColor = 'White';
annot.FontSize = 11;
saveas(figure(1),'gold\OptimalMDot.fig');

table_water = [m',T_prod(1,:)',W_net_IP_MW(1,:)', SpCC_greenfield(1,:)', specific_w(1,:)'];
table_co2 = [m',T_prod(2,:)' ,W_net_IP_MW(2,:)', SpCC_greenfield(2,:)', specific_w(2,:)'];
table = [table_water, table_co2]
table = array2table(table)
writetable(table,'gold\Figure1_Data.xlsx')