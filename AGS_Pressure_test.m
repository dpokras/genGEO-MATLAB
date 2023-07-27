clear;
clc;
close all;

%% Define System Parameters
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
params.wellFieldType = 'Doublet'; 
params.time_years = 30;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.system = 'Conduction4';
params.orcFluid = 'R245fa';
params.res_length = 5000; % in meters
params.depth = 3500; % in meters
params.fluid = ["CO2"];
params.well_radius = 0.25;
%% Determine Optimum m-dot for fluid
result = total_analytic_system_optmdot(params);
params.m_dot_IP = result.m_dot_IP;
temp_at_bottom = params.T_surface_rock_C + params.dT_dz * params.depth;

result = total_analytic_system_co2(params);
z_IW_Water = result.injection.State(:,1)*(-1);
z_LW_Water = (result.reservoir.State(:,1)- params.depth)*(-1);
z_PWu_Water = result.productionUpper.State(:,1);
z_PWu_Water = flip(z_PWu_Water)
T_IW_Water = result.injection.State(:,3);
T_LW_Water = result.reservoir.State(:,3);
T_PWu_Water = result.productionUpper.State(:,3)
P_IW_Water = result.injection.State(:,2);
P_LW_Water = result.reservoir.State(:,2);
P_PWu_Water = result.productionUpper.State(:,2);



%% Plot Figure
%% Plot Pressure Profiles
subplot(1,2,1)
w_P_inj = plot(P_IW_Water,z_IW_Water,'--b');
hold on
w_P_res = plot(P_LW_Water,z_LW_Water,'--b');
w_P_prod = plot(P_PWu_Water,z_PWu_Water,'--b');
title('Pressure Profile for Water')
xlabel('Pressure [MPa]', 'FontSize',11)  % x-axis label
ylabel('Resevoir Depth [km]', 'FontSize',11) % y-axis label
set(gca, 'YDir','reverse')
set(gca,'FontSize',11)
%     % Plot labels for P-profile
%     plot(P_IW_Water(1), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','b', 'MarkerFaceColor', 'w'); grid on; hold on
%     strValues = strtrim(cellstr(num2str(1)));
%     text(P_IW_Water(1),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'b' );
%     plot(P_LW_Water(1), params.depth/1000, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','b', 'MarkerFaceColor', 'w'); grid on; hold on
%     strValues = strtrim(cellstr(num2str(2)));
%     text(P_LW_Water(1),params.depth/1000,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'b');
%     plot(P_LW_Water(end), params.depth/1000, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','b', 'MarkerFaceColor', 'w'); grid on; hold on
%     strValues = strtrim(cellstr(num2str(3)));
%     text(P_LW_Water(end),params.depth/1000,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'b');
%     plot(P_PW_Water(end), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','b', 'MarkerFaceColor', 'w'); grid on; hold on
%     strValues = strtrim(cellstr(num2str(4)));
%     text(P_PW_Water(end),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'b');

%% Plot Temperature Profile
subplot(1,2,2)
w_T_inj = plot(T_IW_Water,z_IW_Water,'-b');
hold on
w_T_res = plot(T_LW_Water,z_LW_Water,'-b');
plot(T_PWu_Water,z_PWu_Water,'-b');
title('Temperature Profile for Water')
xlabel('Temperature [\circ{}C]', 'FontSize',11)  % x-axis label
ylabel('Resevoir Depth [km]', 'FontSize',11) % y-axis label
set(gca, 'YDir','reverse')
set(gca,'FontSize',11)

%     % Plot Labels for T profile
%     plot(T_IW_Water(1), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','b', 'MarkerFaceColor', 'w'); grid on; hold on
%     strValues = strtrim(cellstr(num2str(1)));
%     text(T_IW_Water(1),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'b' );
%     plot(T_LW_Water(1), params.depth/1000, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','b', 'MarkerFaceColor', 'w'); grid on; hold on
%     strValues = strtrim(cellstr(num2str(2)));
%     text(T_LW_Water(1),params.depth/1000,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'b');
%     plot(T_LW_Water(end), params.depth/1000, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','b', 'MarkerFaceColor', 'w'); grid on; hold on
%     strValues = strtrim(cellstr(num2str(3)));
%     text(T_LW_Water(end),params.depth/1000,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'b');
%     plot(T_PW_Water(end), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','b', 'MarkerFaceColor', 'w'); grid on; hold on
%     strValues = strtrim(cellstr(num2str(4)));
%     text(T_PW_Water(end),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'b');
%     
saveas(figure(1),'images\Pressure_Profiles_test.fig');