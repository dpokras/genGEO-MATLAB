clear;
close all;
clc;
% Purpose of sheet is to form a figure out of the .xlsx generated by
% AGS_Num_Laterals.m
params = SimulationParameters;
params.dT_dz = 0.035;
%% Financial Aspects
F_OM = 0.045;
d = 0.096;
n = 25;
CF = 0.85;
CRF = (d*(1+d)^n)/((1+d)^n-1);
RF = (CRF + F_OM) / (CF * 8760) * 1e6
%% Extract Depth and Corresponding Optimal Reservoir Lengths
Conduction_1_data = xlsread('data\Num_RWs_Conduction1.xlsx');
Conduction_2_data = xlsread('data\Num_RWs_Conduction2.xlsx');
Conduction_4_data = xlsread('data\Num_RWs_Conduction4.xlsx');
Conduction_8_data = xlsread('data\Num_RWs_Conduction8.xlsx');
% Data is in form: [Depth, optimum Res-length,optimum m-dot,Power, Specific Power, Prod_T, LCOE]
depth = Conduction_1_data(:,1)/1000;
temp_at_depth = params.dT_dz * depth * 1000 + params.T_surface_rock_C;
%% Extract Data for Base Case
base_case = xlsread('data\Base Case Results.xlsx');
base_case(1) = base_case(1)/1000;
%% Select Variables in Plot
sys = ["Conduction1","Conduction2","Conduction4","Conduction8"];
%% Plot Figures
figure(1)
axis tight
% Optimal Reservoir Length Subplot
subplot(3,2,3)
con1 = plot(NaN,NaN,':k','linewidth',1.4);
hold on
con2 = plot(NaN,NaN,'-.k','linewidth',1.4);
con3 = plot(NaN,NaN,'--k','linewidth',1.4);
con4 = plot(NaN,NaN,'-k','linewidth',1.4);
space = plot(NaN,NaN,'-w','linewidth',1.4);
grad = plot(NaN,NaN,'-','Color','#a8a8a8','linewidth',1.4);
star = plot(NaN,NaN,'pentagram','MarkerFaceColor','#000000','MarkerEdgeColor','#000000');
% opt res length
plot(Conduction_1_data(:,2)/1000,depth,':k','linewidth',1.4);
plot(Conduction_2_data(:,2)/1000,depth,'-.k','linewidth',1.4);
plot(Conduction_4_data(:,2)/1000,depth,'--k','linewidth',1.4);
optreslength = plot(Conduction_8_data(:,2)/1000,depth,'-k','linewidth',1.4);
plot(base_case(2)/1000,base_case(1),'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black');
xlabel('Optimal Lateral Well Length [km]', 'FontSize',10, 'HorizontalAlignment','right'); % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
ax1 = gca;
set(ax1, 'Position', [0.08 0.425 0.4 0.23],'YDir','reverse','Ylim', [min(depth),max(depth)],'Xlim', [0,35],'fontsize',10.5,'xColor','k','xtick',0:5:15,'gridcolor','k');
grid on
ax2 = axes('Position',[0.08 0.425 0.4 0.23],...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
set(ax2,'FontSize',11)
hold(ax2,'on')
% opt m-dot
plot(Conduction_1_data(:,3),depth,':','Color','#4e5fa2','Parent',ax2)
plot(Conduction_2_data(:,3),depth,'-.','Color','#4e5fa2','Parent',ax2)
plot(Conduction_4_data(:,3),depth,'--','Color','#4e5fa2','Parent',ax2)
optmdot = plot(Conduction_8_data(:,3),depth,'-','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(base_case(3),base_case(1),'pentagram','MarkerFaceColor','#4e5fa2','MarkerEdgeColor','#4e5fa2','Parent',ax2);
set(ax2, 'YDir','reverse','Ylim', [min(depth),max(depth)],'Xlim', [-120,160],'fontsize',11,'xColor','#4e5fa2','xtick',0:40:160,'gridcolor','k');
set(ax2,'YTickLabel',[]);
ylim([1,8])

xlabel( 'Optimal Mass Flowrate [kg s^{-1}]', 'FontSize',11,'HorizontalAlignment','left'); % x-axis label
% legend
lgd = legend([con1 con2 con3 con4 star grad],"1 Lateral", "2 Laterals", "4 Laterals", "8 Laterals", "Base Case","Geothermal Gradient",'numcolumns',3);
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 11;
lgd.FontSize = 11;
lgd.Position = [0.15 0.3 0.1 0.2]; %[110 500 150 90];
grid on

% Power Subplot
subplot(3,2,5)
plot(Conduction_1_data(:,4)/2,depth,':k','linewidth',1.4);
hold on
plot(Conduction_2_data(:,4)/2,depth,'-.k','linewidth',1.4);
plot(Conduction_4_data(:,4)/2,depth,'--k','linewidth',1.4);
plot(Conduction_8_data(:,4)/2,depth,'-k','linewidth',1.4);
plot(base_case(4),base_case(1),'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black')
set(gca, 'Position', [0.08 0.07 0.4 0.23],'YDir','reverse')
set(gca,'FontSize',11)
xlabel('Electric Power [MW_{e}]', 'FontSize',11) % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
grid on
ylim([1,8])
% Production Temperature Subplot
subplot(3,2,6)
plot(Conduction_1_data(:,6),depth,':k');
hold on
plot(Conduction_2_data(:,6),depth,'-.k');
plot(Conduction_4_data(:,6),depth,'--k');
plot(Conduction_8_data(:,6),depth,'-k');
plot(temp_at_depth, depth, '-','Color','#a8a8a8','linewidth',1.4);
plot(base_case(6),base_case(1),'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black')
set(gca, 'Position',[0.55 0.07 0.4 0.23],'YDir','reverse')
set(gca,'FontSize',11)
annotation('textbox', [0.75, 0.15, 0.1, 0.1], 'String', "\nabla T = 35 \circ{}C km^{-1}",'Color','#8c8c8c','EdgeColor','#8c8c8c')
xlabel('Production Temperature [\circ{}C]', 'FontSize',11) % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
grid on
ylim([1,8])
% SpCC Subplot and Specific Power
subplot(3,2,4)
plot(Conduction_1_data(:,5),depth,':k','linewidth',1.4);
hold on
plot(Conduction_2_data(:,5),depth,'-.k','linewidth',1.4);
plot(Conduction_4_data(:,5),depth,'--k','linewidth',1.4);
plot(Conduction_8_data(:,5),depth,'-k','linewidth',1.4);
plot(base_case(5),base_case(1),'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black')
set(gca, 'Position', [0.55 0.425 0.4 0.23],'YDir','reverse','Ylim', [1,8],'Xlim', [0,80],'fontsize',11,'xColor','black','gridcolor','k');
set(gca,'FontSize',11)
xlabel('Specific Electric Power [W_{e} m^{-1}]', 'FontSize',11) % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
ax1 = gca;
ax1_pos = ax1.Position;
% spcc
ax2 = axes('Position',[0.55 0.425 0.4 0.23],...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
set(ax2,'FontSize',11)
hold(ax2,'on')
plot(Conduction_1_data(:,7)/RF,depth,':','Color','#4e5fa2','Parent',ax2);
plot(Conduction_2_data(:,7)/RF,depth,'-.','Color','#4e5fa2','Parent',ax2);
plot(Conduction_4_data(:,7)/RF,depth,'--','Color','#4e5fa2','Parent',ax2);
plot(Conduction_8_data(:,7)/RF,depth,'-','Color','#4e5fa2','Parent',ax2);
plot(base_case(7),base_case(1),'pentagram','MarkerFaceColor','#4e5fa2','MarkerEdgeColor','black','Parent',ax2)
set(ax2, 'YDir','reverse','Ylim', [1,8],'Xlim', [0,4000],'fontsize',11,'xColor','#4e5fa2','gridcolor','k');
set(ax2,'YTickLabel',[]);
xlabel('SpCC [$ W_{e}^{-1}]', 'Fontsize',11) % x-axis label
grid on

caption = {'{\bfProperties}'...
    ,'Massflow to Minimize Cost'...
    ,'Lateral Well Length to Minimize Cost'...
    ,'CO_{2} Working Fluid'...
    ,'35 �C km^{-1} Geothermal Gradient'...
    ,'15 �C Surface Temp'...
    ,'0.5 m Well Diameter'...
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




saveas(figure(1),'images\Num_Laterals.fig');