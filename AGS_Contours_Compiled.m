clear;
close all;
clc;

%% Financial Aspects
F_OM = 0.045;
d = 0.096;
n = 25;
CF = 0.85;
CRF = (d*(1+d)^n)/((1+d)^n-1);
RF = (CRF + F_OM) / (CF * 8760) * 1e6;
%% Read Data from .xlsx Files
% LCOE Data
CO2_LCOE_Baseline = xlsread('data\Contour_CO2_LCOE_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
CO2_LCOE_Ideal = xlsread('data\Contour_CO2_LCOE_Conduction4_dTdz35_r25_Ideal.xlsx');
d_co2_baseline = CO2_LCOE_Baseline(:,1)/1000;
d_co2_ideal = CO2_LCOE_Ideal(:,1)/1000;
r_co2_baseline = CO2_LCOE_Baseline(:,2)/1000;
r_co2_ideal = CO2_LCOE_Ideal(:,2)/1000;
CO2_SpCC_Baseline = CO2_LCOE_Baseline(:,3)/RF;
CO2_SpCC_Ideal = CO2_LCOE_Ideal(:,3)/RF;
% Power Data
CO2_Power = xlsread('data\Contour_CO2_Power_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
CO2_Power = CO2_Power(:,3)/2;
% Specific Power Data
CO2_Specific_Power = xlsread('data\Contour_CO2_Specific_Power_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
CO2_Specific_Power = CO2_Specific_Power(:,3);
% Production Temperature Data
CO2_Temp = xlsread('data\Contour_CO2_Production_Temp_Conduction4_dTdz35_r25_Baseline-Feb2021.xlsx');
CO2_Temp = CO2_Temp(:,3);
% Optimum Reservoir Length Data
CO2_Res_Length_Baseline = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz35_radius0.25.xlsx');
CO2_Res_Length_Ideal = xlsread('data\Optimum_Res_Length_CO2_Conduction4_dTdz35_radius25_ideal.xlsx');
opt_res_length_co2_baseline = CO2_Res_Length_Baseline(:,2)/1000;
opt_res_length_co2_ideal = CO2_Res_Length_Ideal(:,2)/1000;
depth_co2_baseline = CO2_Res_Length_Baseline(:,1)/1000;
depth_co2_ideal = CO2_Res_Length_Ideal(:,1)/1000;
%% create gidspace for contour plot
[X,Y] = meshgrid(linspace(min(r_co2_baseline),max(r_co2_baseline),2000), linspace(min(d_co2_baseline),max(d_co2_baseline),2000));
%% Plot LCOE Contour Data
figure(1)

contourLevels = [50, 100, 250, 500, 2500];
% LCOE Plot Baseline
[M1,c1] =contour(X,Y,griddata(r_co2_baseline,d_co2_baseline,CO2_SpCC_Baseline',X,Y));
c1.LevelList = contourLevels;
c1.LineColor = '#000000';
c1.LineStyle = "-";
c1.LineWidth = 1;
c1.DisplayName = 'Baseline';
clabel(M1,c1,'manual');
hold on
opt_baseline = plot(opt_res_length_co2_baseline,depth_co2_baseline,'-','Color','#4e5fa2','linewidth',1.4);
star = plot(5,3.5,'pentagram','MarkerFaceColor','#000000','MarkerEdgeColor','#000000');
%LCOE Plot Ideal
[M2,c2] =contour(X,Y,griddata(r_co2_ideal,d_co2_ideal,CO2_SpCC_Ideal',X,Y));
c2.LevelList = contourLevels;
c2.LineColor = '#000000';
c2.LineStyle = "--";
c2.LineWidth = 1;
c2.DisplayName = 'Ideal';
clabel(M2,c2,'manual');
title('SpCC [$ W_{e}^{-1}]', 'FontSize',11)
xlabel('Lateral Well Length [km]', 'FontSize',11) % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
opt_ideal = plot(opt_res_length_co2_ideal,depth_co2_ideal,'--','Color','#4e5fa2','linewidth',1.4);
lgd = legend([c1,c2,opt_baseline,opt_ideal, star],'SpCC, Baseline Costs','SpCC, Ideal Costs', 'L_{LW}^*, Baseline Costs','L_{LW}^*, Ideal Costs', 'Base Case','numcolumns',3);
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 11;
lgd.FontSize = 11;
xlim([0,30])
ylim([1,8])
set(gca,'FontSize',11, 'Position',  [0.08 0.425 0.4 0.23],'YDir','reverse')

%% Plot Power Contour Data
subplot(3,2,4)
contourLevels = [100,500,1000,2500,5000];
% Power Plot CO2
[M3,c3] =contour(X,Y,griddata(r_co2_baseline,d_co2_baseline,CO2_Power',X,Y));
c3.LevelList = contourLevels;
c3.LineColor = '#000000';
c3.LineWidth = 1;
c3.DisplayName = 'CO_{2}';
clabel(M3,c3,contourLevels,'Margin',5,'LabelSpacing',400,'FontSize',10);
hold on
plot(5,3.5,'pentagram','MarkerFaceColor','#000000','MarkerEdgeColor','#000000');
title('Electric Power [kW_{e}]', 'FontSize',11)
xlabel('Lateral Well Length [km]', 'FontSize',11) % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
xlim([0,30])
ylim([1,8])
set(gca,'FontSize',11, 'Position',  [0.55 0.425 0.4 0.23],'YDir','reverse')

%% Plot Specific Power Contour Data
subplot(3,2,5)
contourLevels = [5,10,20,40,60,80,100,120];
% Specific Power Plot CO2
[M5,c5] =contour(X,Y,griddata(r_co2_baseline,d_co2_baseline,CO2_Specific_Power',X,Y));
c5.LevelList = contourLevels;
c5.LineColor = '#000000';
c5.LineWidth = 1;
c5.DisplayName = 'CO_{2}';
clabel(M5,c5,contourLevels,'Margin',5,'LabelSpacing',400,'FontSize',10);
hold on
plot(5,3.5,'pentagram','MarkerFaceColor','#000000','MarkerEdgeColor','#000000');
title('Specific Electric Power [W_{e} m^{-1}]', 'FontSize',11)
xlabel('Lateral Well Length [km]', 'FontSize',11) % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
xlim([0,30])
ylim([1,8])
set(gca,'FontSize',11,'Position', [0.08 0.07 0.4 0.23],'YDir','reverse')
%% Plot Production Temperature Contour Data
subplot(3,2,6)
contourLevels = [20,40,60,80,100,120,140,160];
% Production Temperature Plot CO2
[M7,c7] =contour(X,Y,griddata(r_co2_baseline,d_co2_baseline,CO2_Temp',X,Y));
c7.LevelList = contourLevels;
c7.LineColor = '#000000';
c7.LineWidth = 1;
c7.DisplayName = 'CO_{2}';
clabel(M7,c7,contourLevels,'Margin',5,'LabelSpacing',400,'FontSize',10);
hold on
plot(5,3.5,'pentagram','MarkerFaceColor','#000000','MarkerEdgeColor','#000000');
title('Production Temperature [\circ{}C]', 'FontSize',11)
xlabel('Lateral Well Length [km]', 'FontSize',11) % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
xlim([0,30])
ylim([1,8])
set(gca,'FontSize',11,'Position', [0.55 0.07 0.4 0.23],'YDir','reverse')

%% Specification Legend
caption = {'{\bfProperties}'...
    ,'Massflow to Minimize Cost'...
    ,'CO_{2} Working Fluid'...
    ,'35 �C km^{-1} Geothermal Gradient'...
    ,'0.5 m Well Diameter'...
    ,'15 �C Surface Temp'...
    ,'Four Horizontal Laterals'...
    ,strcat(['Year-' num2str(30) ' Depletion Values'])...
    ,'Malek|Adams|Rossi|Schiegg|Saar (2021)'};
annot = annotation('textbox','String',caption);
annot.Position = [0.085 0.7 0.3 0.3];
annot.FitBoxToText = 'on';
annot.BackgroundColor = 'White';
annot.FontSize = 11;

annota = annotation('textbox',[0.08 0.68 0.05 0.05],'String',"(a)",'EdgeColor','none','FontSize',14,'fontweight', 'bold');
annotb = annotation('textbox',[0.55 0.68 0.05 0.05],'String',"(b)",'EdgeColor','none','FontSize',14,'fontweight', 'bold');
annotc = annotation('textbox',[0.08 0.31 0.05 0.05],'String',"(c)",'EdgeColor','none','FontSize',14,'fontweight', 'bold');
annotd = annotation('textbox',[0.55 0.31 0.05 0.05],'String',"(d)",'EdgeColor','none','FontSize',14,'fontweight', 'bold');

saveas(figure(1),'images\Contours_Conduction4_Compiled_dTdz35_r25_BaselineIdeal.fig');

