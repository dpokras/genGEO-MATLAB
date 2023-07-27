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

%% Well Radius vs Well Radius Ratio
% radii = linspace(0.2,1.0,5); %m
% radius_ratios = linspace(0.25,1.0,4);
% A = zeros(size(radii,2),17);
% S = struct;
% 
% Z = zeros(size(radii,2),size(radius_ratios,2));
% SpCC = zeros(size(radii,2),size(radius_ratios,2));
% M = zeros(size(radii,2),size(radius_ratios,2));
% C_Gr = zeros(size(radii,2),size(radius_ratios,2));
% C_Surf = zeros(size(radii,2),size(radius_ratios,2));
% C_Sub = zeros(size(radii,2),size(radius_ratios,2));
% W_net = zeros(size(radii,2),size(radius_ratios,2));
% for i = 1:size(radii,2)
%     
%     filename = strcat('dpo_opt_config1_wellradius_',num2str(radii(i)),'_22_06_23.csv');
%     T = readtable(filename);
%     S(i).T = T;
%     SpCC(i,:) = T.SpCC_W_Gr;
%     M(i,:) = T.m_dot;
%     C_Gr(i,:) = T.C_Gr;
%     C_Surf(i,:) = T.C_SurfacePlant;
%     C_Sub(i,:) = T.C_Subsurface_greenfield;
%     W_net(i,:) = T.W_net;
% end
% 
% figure(1);
% y = radii;
% x = radius_ratios;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,SpCC);
% yticks(y);
% xticks(x);
% ylabel('Well Radius (m)');
% xlabel('Well Radius Ratio (-)');
% cb = colorbar;
% cb.Label.String = 'SpCC Greenfield ($/W)';
% figure(2);
% y = radii;
% x = radius_ratios;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,W_net/1e3);
% yticks(y);
% xticks(x);
% ylabel('Well Radius (m)');
% xlabel('Well Radius Ratio (-)');
% colormap hot
% cb = colorbar;
% cb.Label.String = 'Net Electric Power Generation (MW)';
% 
% figure(3);
% y = radii;
% x = radius_ratios;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,C_Gr/1e6);
% yticks(y);
% xticks(x);
% ylabel('Well Radius (m)');
% xlabel('Well Radius Ratio (-)');
% colormap cool
% cb = colorbar;
% cb.Label.String = 'Total Capital Investment Greenfield (M$)';

% %% Reservoir Length vs Well Radius
% radii = linspace(0.2,1.0,5); %m
% length = 5000:1000:12000;
% A = zeros(size(radii,2),17);
% S = struct;
% 
% Z = zeros(size(length,2),size(radii,2));
% SpCC = zeros(size(length,2),size(radii,2));
% M = zeros(size(length,2),size(radii,2));
% C_Gr = zeros(size(length,2),size(radii,2));
% C_Surf = zeros(size(length,2),size(radii,2));
% C_Sub = zeros(size(length,2),size(radii,2));
% W_net = zeros(size(length,2),size(radii,2));
% for i = 1:size(length,2)
%     
%     %this speeds up the optimisation by bringing the starting value closer
%     %to the optimal    
%     filename = strcat('dpo_opt_config1_res',num2str(length(i)),'_vs_radius_24_06_23','.csv');
%     T = readtable(filename);
%     S(i).T = T;
%     SpCC(i,:) = T.SpCC_W_Gr;
%     M(i,:) = T.m_dot;
%     C_Gr(i,:) = T.C_Gr;
%     C_Surf(i,:) = T.C_SurfacePlant;
%     C_Sub(i,:) = T.C_Subsurface_greenfield;
%     W_net(i,:) = T.W_net;
% end
% 
% figure(1);
% y = length/1e3;
% x = radii;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,SpCC);
% yticks(y);
% xticks(x);
% ylabel('Reservoir Length [km]');
% xlabel('Well Radius [m]');
% cb = colorbar;
% cb.Label.String = 'SpCC Greenfield ($/W)';
% 
% figure(2);
% y = length/1e3;
% x = radii;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,W_net/1e3);
% yticks(y);
% xticks(x);
% ylabel('Reservoir Length [km]');
% xlabel('Well Radius [m]');
% colormap hot
% cb = colorbar;
% cb.Label.String = 'Net Electric Power Generation (MW)';
% 
% figure(3);
% y = length/1e3;
% x = radii;
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,C_Gr/1e6);
% yticks(y);
% xticks(x);
% ylabel('Reservoir Length [km]');
% xlabel('Well Radius [m]');
% colormap cool
% cb = colorbar;
% cb.Label.String = 'Total Capital Investment Greenfield (M$)';

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

%% S_ratio vs res_length vs n_streams
ratios = linspace(0, 1, 6);
streams = 5:10;
length = 5000:1000:10000;
S = struct;
j = 5;
k = 10;

SpCC = zeros(size(streams,2),size(length,2),size(ratios,2));
M = zeros(size(streams,2),size(length,2),size(ratios,2));
C_Gr = zeros(size(streams,2),size(length,2),size(ratios,2));
C_Surf = zeros(size(streams,2),size(length,2),size(ratios,2));
C_Sub = zeros(size(streams,2),size(length,2),size(ratios,2));
W_net = zeros(size(streams,2),size(length,2),size(ratios,2));

y = streams;
x = length/1e3;
[X,Y] = meshgrid(x,y);
figure(1);
hold on

for r = 1:size(ratios,2)
    for n = 1:size(streams, 2)
        filename = strcat('dpo_opt_config3_S_ratio_',num2str(ratios(r)),'_vs_n_streams',num2str(streams(n)),'_vs_res_length_27_06_23.csv');
    
        T = readtable(filename);
        S(n).T = T;
        SpCC(:,n,r) = T.SpCC_dH_Gr(1:size(streams,2),:);
        M(:,n,r) = T.m_dot(1:size(streams,2),:);
        C_Gr(:,n,r) = T.C_Subsurface_greenfield(1:size(streams,2),:);
        C_Surf(:,n,r) = T.C_SurfacePlant(1:size(streams,2),:);
        C_Sub(:,n,r) = T.C_Subsurface_greenfield(1:size(streams,2),:);
        W_net(:,n,r) = T.dH_sold_avg (1:size(streams,2),:);

    end 
    [~,h] = contourf(X,Y,SpCC(:,:,r),15, FaceAlpha=0.6, EdgeColor='black', EdgeAlpha=0.6);   % plot contour at the bottom
    h.ContourZLevel = ratios(r);
    
end
hold off
set(gca, 'dataaspectratio', [1 1 0.25], 'projection', 'perspective', 'box', 'on')
yticks(y);
xticks(x);
ylabel('Number of Subsurface Streams');
xlabel('Reservoir Length per Stream [km]');
zlabel('S Ratio')
colormap jet
cb = colorbar;
cb.Label.String = 'SpCC Greenfield ($/W)';
view(3);

h = rotate3d;
set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
set(gcf, 'ResizeFcn', @align_axislabel)
align_axislabel([], gca)
axislabel_translation_slider;
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

figure(5);
x = ratios;
color_map = jet(size(length,2));
hold on

for l = 1:size(length,2)
    y_values = zeros(size(streams, 2), numel(x));
    
    for n = 1:size(streams,2)
        y_values(n,:) = reshape(SpCC(n,l,:), 1, numel(x));
    end
    
    y_min = min(y_values);
    y_max = max(y_values);

    fill([x, fliplr(x)], [y_max, fliplr(y_min)], color_map(l,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', strcat('l = ', num2str(length(l)/1e3),'km'));
end
xlabel('S Ratio (-)')
ylabel('SpCC Greenfield ($/W)')
legend('Location','northwest')


