clc;
clear;
close;
% %% res_ength vs n_streams 
% length = 5000:1000:12000;
% S = struct;
% j = 5;
% k = 10;
% 
% SpCC = zeros(size(j:k,2),size(length,2));
% M = zeros(size(j:k,2),size(length,2));
% C_Gr = zeros(size(j:k,2),size(length,2));
% C_Surf = zeros(size(j:k,2),size(length,2));
% C_Sub = zeros(size(j:k,2),size(length,2));
% W_net = zeros(size(j:k,2),size(length,2));
% 
% for l = 1:size(length,2)
%     filename = strcat('dpo_opt_config1_res',num2str(length(l)),'_22_06_23.csv');
% 
%     T = readtable(filename);
%     S(l).T = T;
%     SpCC(:,l) = T.SpCC_W_Gr(1:size(j:k,2),:);
%     M(:,l) = T.m_dot(1:size(j:k,2),:);
%     C_Gr(:,l) = T.C_Subsurface_greenfield(1:size(j:k,2),:);
%     C_Surf(:,l) = T.C_SurfacePlant(1:size(j:k,2),:);
%     C_Sub(:,l) = T.C_Subsurface_greenfield(1:size(j:k,2),:);
%     W_net(:,l) = T.W_net(1:size(j:k,2),:);
% 
% end
% 
% figure(1);
% y = j:k;
% x = length/1e3;
% [X,Y] = meshgrid(x,y);
% surf(X,Y,SpCC);
% yticks(y);
% xticks(x);
% ylabel('Number of Subsurface Streams');
% xlabel('Reservoir Length per Stream [km]');
% cb = colorbar;
% cb.Label.String = 'SpCC Greenfield ($/W)';
% 
% figure(2);
% y = j:k;
% x = length/1e3;
% hold on
% for i = flip(1:size(j:k,2))
%     plot(x,M(i,:), 'DisplayName',strcat('n = ', num2str(y(i))))
% end
% ylabel('Optimal Mass Flowrate (kg/s)');
% xlabel('Reservoir Length per Stream [km]');
% legend show
% hold off
% 
% figure(3);
% y = j:k;
% x = length/1e3;
% hold on
% for i = flip(1:size(j:k,2))
%     plot(x,C_Surf(i,:)/1e6, 'DisplayName',strcat('n = ', num2str(y(i))))
% end
% ylabel('Total Capital Investment (Surface) (M$)');
% xlabel('Reservoir Length per Stream [km]');
% legend show
% hold off
% 
% figure(4);
% y = j:k;
% x = length/1e3;
% hold on
% for i = flip(1:size(j:k,2))
%     plot(x,C_Sub(i,:)/1e6, 'DisplayName',strcat('n = ', num2str(y(i))))
% end
% ylabel('Total Capital Investment (Sub-surface) (M$)');
% xlabel('Reservoir Length per Stream [km]');
% legend show
% 
% figure(5);
% y = j:k;
% x = length/1e3;
% hold on
% for i = flip(1:size(j:k,2))
%     plot(x,W_net(i,:)/1e3, 'DisplayName',strcat('n = ', num2str(y(i))))
% end
% ylabel('Net Electric Power Generation (kW)');
% xlabel('Reservoir Length per Stream [km]');
% legend show
% hold off
% 
% figure(6);
% y = j:k;
% x = length/1e3;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,W_net/1e3,15);
% yticks(y);
% xticks(x);
% ylabel('Number of Subsurface Streams');
% xlabel('Reservoir Length per Stream [km]');
% colormap hot
% cb = colorbar;
% cb.Label.String = 'Net Electric Power Generation (MW)';
% 
% figure(7);
% y = j:k;
% x = length/1e3;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,C_Gr/1e6,15);
% yticks(y);
% xticks(x);
% ylabel('Number of Subsurface Streams');
% xlabel('Reservoir Length per Stream [km]');
% colormap cool
% cb = colorbar;
% cb.Label.String = 'Total Capital Investment Greenfield (M$)';
% 
% figure(8);
% y = j:k;
% x = length/1e3;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,C_Sub/1e6,15);
% yticks(y);
% xticks(x);
% ylabel('Number of Subsurface Streams');
% xlabel('Reservoir Length per Stream [km]');
% colormap cool
% cb = colorbar;
% cb.Label.String = 'Total Capital Investment Greenfield (Sub-surface) (M$)';
% 
% figure(9);
% y = j:k;
% x = length/1e3;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,C_Surf/1e6,15);
% yticks(y);
% xticks(x);
% ylabel('Number of Subsurface Streams');
% xlabel('Reservoir Length per Stream [km]');
% colormap cool
% cb = colorbar;
% cb.Label.String = 'Total Capital Investment (Surface) (M$)';

% % %% Well Radius vs Well Radius Ratio
% % radii = linspace(0.2,1.0,5); %m
% % radius_ratios = linspace(0.25,1.0,4);
% % A = zeros(size(radii,2),17);
% % S = struct;
% % 
% % Z = zeros(size(radii,2),size(radius_ratios,2));
% % SpCC = zeros(size(radii,2),size(radius_ratios,2));
% % M = zeros(size(radii,2),size(radius_ratios,2));
% % C_Gr = zeros(size(radii,2),size(radius_ratios,2));
% % C_Surf = zeros(size(radii,2),size(radius_ratios,2));
% % C_Sub = zeros(size(radii,2),size(radius_ratios,2));
% % W_net = zeros(size(radii,2),size(radius_ratios,2));
% % for i = 1:size(radii,2)
% %     
% %     filename = strcat('dpo_opt_config1_wellradius_',num2str(radii(i)),'_22_06_23.csv');
% %     T = readtable(filename);
% %     S(i).T = T;
% %     SpCC(i,:) = T.SpCC_W_Gr;
% %     M(i,:) = T.m_dot;
% %     C_Gr(i,:) = T.C_Gr;
% %     C_Surf(i,:) = T.C_SurfacePlant;
% %     C_Sub(i,:) = T.C_Subsurface_greenfield;
% %     W_net(i,:) = T.W_net;
% % end
% % 
% % figure(1);
% % y = radii;
% % x = radius_ratios;
% % [X,Y] = meshgrid(x,y);
% % contourf(X,Y,SpCC);
% % yticks(y);
% % xticks(x);
% % ylabel('Well Radius (m)');
% % xlabel('Well Radius Ratio (-)');
% % cb = colorbar;
% % cb.Label.String = 'SpCC Greenfield ($/W)';
% % figure(2);
% % y = radii;
% % x = radius_ratios;
% % [X,Y] = meshgrid(x,y);
% % contourf(X,Y,W_net/1e3);
% % yticks(y);
% % xticks(x);
% % ylabel('Well Radius (m)');
% % xlabel('Well Radius Ratio (-)');
% % colormap hot
% % cb = colorbar;
% % cb.Label.String = 'Net Electric Power Generation (MW)';
% % 
% % figure(3);
% % y = radii;
% % x = radius_ratios;
% % [X,Y] = meshgrid(x,y);
% % contourf(X,Y,C_Gr/1e6);
% % yticks(y);
% % xticks(x);
% % ylabel('Well Radius (m)');
% % xlabel('Well Radius Ratio (-)');
% % colormap cool
% % cb = colorbar;
% % cb.Label.String = 'Total Capital Investment Greenfield (M$)';
% % 
% % %% Reservoir Length vs Well Radius
% % radii = linspace(0.2,1.0,5); %m
% % length = 5000:1000:12000;
% % A = zeros(size(radii,2),17);
% % S = struct;
% % 
% % Z = zeros(size(length,2),size(radii,2));
% % SpCC = zeros(size(length,2),size(radii,2));
% % M = zeros(size(length,2),size(radii,2));
% % C_Gr = zeros(size(length,2),size(radii,2));
% % C_Surf = zeros(size(length,2),size(radii,2));
% % C_Sub = zeros(size(length,2),size(radii,2));
% % W_net = zeros(size(length,2),size(radii,2));
% % for i = 1:size(length,2)
% %     
% %     %this speeds up the optimisation by bringing the starting value closer
% %     %to the optimal    
% %     filename = strcat('dpo_opt_config1_res',num2str(length(i)),'_vs_radius_24_06_23','.csv');
% %     T = readtable(filename);
% %     S(i).T = T;
% %     SpCC(i,:) = T.SpCC_W_Gr;
% %     M(i,:) = T.m_dot;
% %     C_Gr(i,:) = T.C_Gr;
% %     C_Surf(i,:) = T.C_SurfacePlant;
% %     C_Sub(i,:) = T.C_Subsurface_greenfield;
% %     W_net(i,:) = T.W_net;
% % end
% % 
% % figure(1);
% % y = length/1e3;
% % x = radii;
% % [X,Y] = meshgrid(x,y);
% % contourf(X,Y,SpCC);
% % yticks(y);
% % xticks(x);
% % ylabel('Reservoir Length [km]');
% % xlabel('Well Radius [m]');
% % cb = colorbar;
% % cb.Label.String = 'SpCC Greenfield ($/W)';
% % 
% % figure(2);
% % y = length/1e3;
% % x = radii;
% % [X,Y] = meshgrid(x,y);
% % contourf(X,Y,W_net/1e3);
% % yticks(y);
% % xticks(x);
% % ylabel('Reservoir Length [km]');
% % xlabel('Well Radius [m]');
% % colormap hot
% % cb = colorbar;
% % cb.Label.String = 'Net Electric Power Generation (MW)';
% % 
% % figure(3);
% % y = length/1e3;
% % x = radii;
% % [X,Y] = meshgrid(x,y);
% % contourf(X,Y,C_Gr/1e6);
% % yticks(y);
% % xticks(x);
% % ylabel('Reservoir Length [km]');
% % xlabel('Well Radius [m]');
% % colormap cool
% % cb = colorbar;
% % cb.Label.String = 'Total Capital Investment Greenfield (M$)';
% 
% %% S_ratio vs res_length
% ratios = linspace(0,1.0,10);
% n = 5:10;
% A = zeros(size(ratios,2),17);
% S = struct;
% 
% Z = zeros(size(n,2),size(ratios,2));
% SpCC_dH = zeros(size(n,2),size(ratios,2));
% M = zeros(size(n,2),size(ratios,2));
% C_Gr = zeros(size(n,2),size(ratios,2));
% C_Surf = zeros(size(n,2),size(ratios,2));
% C_Sub = zeros(size(n,2),size(ratios,2));
% dH_sold = zeros(size(n,2),size(ratios,2));
% 
% 
% for i = 1:size(ratios,2)
%     filename = strcat('dpo_opt_config3_S_ratio_1_24_06_23_.csv');
%     T = readtable(filename);
%     S(i).T = T;
%     SpCC_dH(i,:) = T.SpCC_dH_Gr;
%     M(i,:) = T.m_dot;
%     C_Gr(i,:) = T.C_Subsurface_greenfield;
%     C_Surf(i,:) = T.C_SurfacePlant;
%     C_Sub(i,:)= T.C_Subsurface_greenfield;
%     dH_sold(i,:) = T.dH_sold_avg;
% end
%     
% figure(1);
% plot(ratios, T.SpCC_dH_Gr,  'k','LineWidth', 1);
% ylabel('SpCC Greenfield ($/W)');
% xlabel('S Ratio (-)');
% 
% figure(2);
% plot(ratios, T.dH_sold_avg,  'ko','LineWidth', 1);
% ylabel('Total Power Sold (MW)');
% xlabel('S Ratio (-)');

% %% Thermal Drawdown Power over time
% 
% filename = 'dpo_opt_config1_temperature_drawdown_26_06_23.csv';
% T = readtable(filename);
%     
% figure(1);
% plot(T.year, T.Power_per_year/1e3,'ko', 'LineWidth',1);
% hold on 
% plot([0,T.year.'], [T.dH_sold_avg(1),T.dH_sold_avg.']/1e3, 'r', 'LineWidth',1);
% hold off
% ylabel('Electical Power Generation [MW]');
% xlabel('time [yr]');
% ylim([0,3])
% 
% figure(2);
% plot(T.year, T.m_dot,'ko', 'LineWidth',1);
% ylabel('Mass Flow Rate (kg/s)');
% xlabel('time [yr]');
% %% S_ratio vs res_length vs n_streams
% ratios = linspace(0, 1, 6);
% streams = 5:10;
% length = 5000:1000:10000;
% S = struct;
% j = 5;
% k = 10;
% 
% SpCC = zeros(size(streams,2),size(length,2),size(ratios,2));
% M = zeros(size(streams,2),size(length,2),size(ratios,2));
% C_Gr = zeros(size(streams,2),size(length,2),size(ratios,2));
% C_Surf = zeros(size(streams,2),size(length,2),size(ratios,2));
% C_Sub = zeros(size(streams,2),size(length,2),size(ratios,2));
% W_net = zeros(size(streams,2),size(length,2),size(ratios,2));
% 
% y = streams;
% x = length/1e3;
% [X,Y] = meshgrid(x,y);
% figure(1);
% hold on
% 
% for r = 1:size(ratios,2)
%     for n = 1:size(streams, 2)
%         filename = strcat('dpo_opt_config3_S_ratio_',num2str(ratios(r)),'_vs_n_streams',num2str(streams(n)),'_vs_res_length_27_06_23.csv');
%     
%         T = readtable(filename);
%         S(n).T = T;
%         SpCC(:,n,r) = T.SpCC_dH_Gr(1:size(streams,2),:);
%         M(:,n,r) = T.m_dot(1:size(streams,2),:);
%         C_Gr(:,n,r) = T.C_Subsurface_greenfield(1:size(streams,2),:);
%         C_Surf(:,n,r) = T.C_SurfacePlant(1:size(streams,2),:);
%         C_Sub(:,n,r) = T.C_Subsurface_greenfield(1:size(streams,2),:);
%         W_net(:,n,r) = T.dH_sold_avg (1:size(streams,2),:);
% 
%     end 
%     [~,h] = contourf(X,Y,SpCC(:,:,r),15, FaceAlpha=0.6, EdgeColor='black', EdgeAlpha=0.6);   % plot contour at the bottom
%     h.ContourZLevel = ratios(r);
%     
% end
% hold off
% set(gca, 'dataaspectratio', [1 1 0.25], 'projection', 'perspective', 'box', 'on')
% yticks(y);
% xticks(x);
% ylabel('Number of Subsurface Streams');
% xlabel('Reservoir Length per Stream [km]');
% zlabel('S Ratio')
% colormap jet
% cb = colorbar;
% cb.Label.String = 'SpCC Greenfield ($/W)';
% view(3);
% 
% h = rotate3d;
% set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
% set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
% set(gcf, 'ResizeFcn', @align_axislabel)
% align_axislabel([], gca)
% axislabel_translation_slider;
% 
% figure(2);
% y = streams;
% x = length/1e3;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,SpCC(:,:,1),15);
% yticks(y);
% xticks(x);
% ylabel('Number of Subsurface Streams');
% xlabel('Reservoir Length per Stream [km]');
% cb = colorbar;
% cb.Label.String = 'SpCC Greenfield ($/W)';
% 
% figure(3);
% y = streams;
% x = length/1e3;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,SpCC(:,:,1),15);
% yticks(y);
% xticks(x);
% ylabel('Number of Subsurface Streams');
% xlabel('Reservoir Length per Stream [km]');
% cb = colorbar;
% cb.Label.String = 'SpCC Greenfield ($/W)';
% 
% 
% figure(5);
% x = ratios;
% color_map = jet(size(length,2));
% hold on
% 
% for l = 1:size(length,2)
%     y_values = zeros(size(streams, 2), numel(x));
%     
%     for n = 1:size(streams,2)
%         y_values(n,:) = reshape(SpCC(n,l,:), 1, numel(x));
%     end
%     
%     y_min = min(y_values);
%     y_max = max(y_values);
% 
%     fill([x, fliplr(x)], [y_max, fliplr(y_min)], color_map(l,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', strcat('l = ', num2str(length(l)/1e3),'km'));
% end
% xlabel('S Ratio (-)')
% ylabel('SpCC Greenfield ($/W)')
% legend('Location','northwest')

% %% S on W_el and W_th and dH
% % Assuming S_ratio, Total_power_avg, Q_net_avg, and W_net_avg are in T
% filename = 'Results/dpo_opt_config3_S_ratio_power_variance_02_08_23.csv';
% 
% T = readtable(filename);
% 
% % Create a new figure
% figure;
% 
% % Plot ΔH, Qth, and Wel vs S_ratio
% plot(T.S_ratio(2:end), T.Total_power_avg(2:end) / 1e6, 'k', 'LineWidth', 2, 'LineStyle', '-'); % ΔH
% hold on;
% plot(T.S_ratio(2:end), T.Q_net_avg(2:end) / 1e6, 'k', 'LineWidth', 2, 'LineStyle', '--'); % Q_th
% plot(T.S_ratio(2:end), T.W_net_avg(2:end) / 1e6, 'k', 'LineWidth', 2, 'LineStyle', ':'); % W_el
% plot(T.S_ratio(1), T.Total_power_avg(1) / 1e6, 'k', 'LineWidth', 2, 'Marker', 'o','MarkerFaceColor','k'); % ΔH
% 
% % Add vertical dotted red line
% line([0.1691 0.1691], ylim, 'Color', 'r', 'LineStyle', '--','LineWidth', 1.6);
% 
% % Label the axes
% xlabel('Mass Flow Split Ratio (-)');
% ylabel('Power Generation (MW)');
% xticks(0:0.1:1);
% 
% % Create a legend and enlarge the font
% lgd = legend('\DeltaH', 'Q_{th}', 'W_{el}', 'Location', 'best');
% lgd.FontSize = 12;
% 
% % Add a grid for better visibility
% grid off;
% 
% % Release the figure
% hold off;
% 
% figure;
% % Plot Mdot vs S_ratio
% plot(T.S_ratio(2:end), T.m_dot(2:end), 'k', 'LineWidth', 2, 'LineStyle', '-'); % ΔH
% hold on;
% plot(T.S_ratio(1), T.m_dot(1), 'k', 'LineWidth', 2, 'Marker', 'o','MarkerFaceColor','k'); % ΔH
% 
% % Add vertical dotted red line
% line([0.1691 0.1691], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.6);
% 
% xlabel('Mass Flow Split Ratio (-)');
% ylabel('Mass Flow Rate (kg/s)');
% xticks(0:0.1:1);
% 
% grid off;
% 
% hold off;
% %% Config Comparison
% % File names for the 4 configurations
% fileNames = {'config1', 'config2', 'config3_S0.5', 'config4'};
% 
% % Configurations names for tick labels
% configs = {'1','2','3 (S = 0.5)','4'};
% 
% % Empty array to store SpCC values
% SpCC = zeros(1, length(fileNames));
% 
% % Loop through the file names
% for i = 1:length(fileNames)
%     % Create the file path
%     filePath = strcat('Results/dpo_', fileNames{i}, '_base_SpCC_30_07_23.csv');
%     
%     % Read the CSV file
%     data = readtable(filePath);
%     
%     % Get the SpCC values from the specific column
%     SpCC(i) = mean(data.SpCC_dH_Gr);
% end
% 
% % Define hatch angles
% hatch_angles = [15, -45, 60, -60]; 
% hatch_density = [15,10,5,5];
% 
% % Create figure
% figure;
% hold on;
% barCenters = linspace(1, length(SpCC), length(SpCC));
% barWidth = 0.8; % This can be adjusted
% 
% hatchOffset = 0.5; % Amount to offset the hatch filling from the x-axis. Adjust this as necessary.
% 
% for i = 1:length(SpCC)
%     % Create patch object for each bar
%     p = patch([barCenters(i)-barWidth/2, barCenters(i)+barWidth/2, barCenters(i)+barWidth/2, barCenters(i)-barWidth/2], ...
%               [hatchOffset, hatchOffset, SpCC(i), SpCC(i)], 'w'); % Offset the hatch filling from the x-axis
%     
%     % Apply hatch pattern
%     hatchfill(p, 'single', hatch_angles(i), hatch_density(i));
% end
% 
% % Set the tick positions and labels
% xticks(barCenters);
% xticklabels(configs);
% 
% % Label axes
% xlabel('Configuration');
% ylabel('SpCC_{\Delta H} Greenfield (2022$/ W)');
% 
% ylim([0,300])
% 
% hold off;

%% Config 1 M dot variance
% % This script reads data from a CSV file, extracts the required data, and
% % then plots SpCC Greenfield vs mass flow rate on the left y-axis, and
% % Electrical Power Generation vs mass flow rate on the right y-axis.
% 
% % Read the CSV file
% resultsTable = readtable('Results/dpo_opt_config1_m_dot_variance_31_07_23.csv');
% 
% % Extract the required data
% m_dot = resultsTable.m_dot;
% SpCC_Gr = resultsTable.SpCC_Gr;
% Power = resultsTable.PowerTotalSoldAvg/1e6;
% 
% % Create a new figure
% figure;
% 
% % Plot SpCC_Gr vs m_dot on the left y-axis
% yyaxis left;
% plot(m_dot, SpCC_Gr, 'k', 'LineWidth', 1.5); % Thinner, black line
% ylabel('SpCC Greenfield (2022$/W)', 'FontSize', 12, 'Color', 'k');
% ylim([0, 1500]); % Set the limits of the y-axis
% yticks(0:250:1500); % Set the intervals of the y-axis
% 
% % Change color of the left y-axis ticks to black
% set(gca,'ycolor','k');
% 
% % Plot Power vs m_dot on the right y-axis
% yyaxis right;
% plot(m_dot, Power, '--k', 'LineWidth', 1.5); % Thinner, black line
% ylim([0, 2]); % Set the limits of the y-axis
% ylabel('Net Electrical Power Generation (MW_{el})', 'FontSize', 12, 'Color', 'k');
% xlabel('Mass Flow Rate (kg/s)', 'FontSize', 12, 'Color', 'k');
% 
% xlim([0, 190]); % Set the limits of the x-axis
% 
% % Change color of the right y-axis ticks to black
% ax = gca;
% ax.YAxis(2).Color = 'k';
% 
% % Add a legend
% legend('SpCC_{Gr}', 'Power', 'FontSize', 12, 'Location', 'best', 'TextColor', 'k');
% grid on;
% 

%% Cost Breakdown
% Load the data
T = readtable('Results/dpo_config3_S05_cost_breakdown_03_08_23.csv');

% Variables specific to Greenfield (Gr)
totalGreenField = [T.C_wells_horizontal, T.C_wells_vertical_Gr,T.C_TBM_plant +  T.C_tCO2, T.C_contingency_Gr, T.C_WC_Gr, T.C_wellfield];

% Variables specific to Brownfield (Br)
totalBrownField = [T.C_wells_horizontal, T.C_wells_vertical_Br,T.C_TBM_plant +  T.C_tCO2, T.C_contingency_Br, T.C_WC_Br, T.C_wellfield];

% Labels for the costs
labels = {'C_{wells, h}','C_{wells,v}','C_{plant}','C_{contingency}', 'C_{WC}', 'C_{wellfield}' };

% Color map
cmap = flipud(gray(7)); % 7 to match the number of categories

value_gr = compose(['%.3f $M'], totalGreenField/1e6);
percents_gr = compose([newline '%.3f%%'], totalGreenField/sum(totalGreenField)*100);

value_br = compose(['%.3f $M'], totalBrownField/1e6);
percents_br = compose([newline '%.3f%%'], totalBrownField/sum(totalBrownField)*100);

% Pie chart for Greenfield scenario
figure;
pie(totalGreenField, strcat(value_gr, percents_gr));
colormap(cmap);
lgd = legend(labels, 'Location','northwest')
lgd.Position = [0.05, 0.75, 0.2, 0.1];

text(-1.8, -1.1, strcat('Total Cost of Capital:' ,{' '},num2str(sum(totalGreenField)/1e6), '$M'));
% Pie chart for Brownfield scenario
figure;
pie(totalBrownField, strcat(value_br, percents_br));
colormap(cmap);
lgd = legend(labels, 'Location','northwest')
lgd.Position = [0.05, 0.75, 0.2, 0.1];

text(-1.8, -1.1, strcat('Total Cost of Capital:' ,{' '},num2str(sum(totalBrownField)/1e6), '$M'));

