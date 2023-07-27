clear;
close all;
clc;

%% Financial Aspects
F_OM = 0.045;
d = 0.096;
n = 25;
CF = 0.85;
CRF = (d*(1+d)^n)/((1+d)^n-1);
RF = (CRF + F_OM) / (CF * 8760) * 1e6
%% Extract Depth and Corresponding Optimal Reservoir Lengths
radius_1_data = xlsread('data\radius_4.xlsx');
radius_2_data = xlsread('data\radius_6.xlsx');
radius_3_data = xlsread('data\radius_9-5_8.xlsx');
radius_4_data = xlsread('data\radius_basecase.xlsx');
% Data is in form: [Depth, optimum Res-length,optimum m-dot,Power, Specific Power, Prod_T, LCOE]
depth = radius_1_data(:,1)/1000;
%% Extract Data for Base Case
base_case = xlsread('data\Base Case Results.xlsx');
base_case(1) = base_case(1)/1000;
%% Plot Figures
figure(1)
% Optimal Reservoir Length and opimal m-dot Subplot
subplot(2,2,3)
%reservoir lengths
axis tight
plot(radius_1_data(:,2)/1000,depth,':k','linewidth',1.4);
hold on
plot(radius_2_data(:,2)/1000,depth,'-.k','linewidth',1.4);
plot(radius_3_data(:,2)/1000,depth,'--k','linewidth',1.4);
plot(radius_4_data(:,2)/1000,depth,'-k','linewidth',1.4);
plot(base_case(2)/1000,base_case(1),'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black');
xlabel('Optimal Lateral Well Length [km]', 'HorizontalAlignment','right'); % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
set(gca,'FontSize',11,'xColor','black')
grid on
ax1 = gca;
set(ax1,  'Position', [0.06 0.13 0.40 0.5], 'YDir','reverse','Ylim', [min(depth),max(depth)], 'Xlim', [0,16],'xtick',0:2:8,'gridcolor','k','fontsize',10.5);
ax2 = axes('Position',[0.06 0.13 0.40 0.5],...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
set(ax2,'FontSize',11)
hold(ax2,'on')
%opt m-dot
plot(radius_1_data(:,3),depth,':','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(radius_2_data(:,3),depth,'-.','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(radius_3_data(:,3),depth,'--','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(radius_4_data(:,3),depth,'-','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(base_case(3),base_case(1),'pentagram','MarkerFaceColor','#4e5fa2','MarkerEdgeColor','#4e5fa2','Parent',ax2);
set(ax2, 'YDir','reverse','Ylim', [min(depth),max(depth)], 'Xlim', [-140,140],'xColor','#4e5fa2','xtick',0:35:140,'gridcolor','k');
set(ax2,'YTickLabel',[]);
ylim([1,8])
xlabel('Optimal Mass Flowrate [kg s^{-1}]', 'FontSize',11, 'HorizontalAlignment','left'); % x-axis label
grid on

% SpCC and Power Subplot
subplot(2,2,[2,4])
con1 = plot(NaN,NaN,':k','linewidth',1.4);
hold on
con2 = plot(NaN,NaN,'-.k','linewidth',1.4);
con3 = plot(NaN,NaN,'--k','linewidth',1.4);
con4 = plot(NaN,NaN,'-k','linewidth',1.4);
space = plot(NaN,NaN,'-w','linewidth',1.4);
optreslength = plot(NaN,NaN,'-m','linewidth',1.4);
optmdot = plot(NaN,NaN,'-b','linewidth',1.4);
star = plot(NaN,NaN,'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black');
%power plots
plot(radius_1_data(:,4)/2,depth,':k','linewidth',1.4);
plot(radius_2_data(:,4)/2,depth,'-.k','linewidth',1.4);
plot(radius_3_data(:,4)/2,depth,'--k','linewidth',1.4);
pow = plot(radius_4_data(:,4)/2,depth,'-k','linewidth',1.4);
plot(base_case(4),base_case(1),'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black');
set(gca, 'YDir','reverse','Ylim', [min(depth),max(depth)])
set(gca,'FontSize',11)
xlabel('Electric Power [MW_{e}]', 'FontSize',11) % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
set(gca,  'Position', [0.57 0.13 0.40 0.5], 'YDir','reverse','Ylim', [min(depth),max(depth)],'xlim',[0,3],'xColor','black','xtick',0:0.5:3,'gridcolor','k');
ax3 = gca;
ax3_pos = ax3.Position;
ax4 = axes('Position',ax3_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
set(ax4,'FontSize',11)
hold(ax4,'on')
% SpCC plots
plot(radius_1_data(:,7)/RF,depth,':','Color','#4e5fa2','Parent',ax4,'linewidth',1.4);
plot(radius_2_data(:,7)/RF,depth,'-.','Color','#4e5fa2','Parent',ax4,'linewidth',1.4);
plot(radius_3_data(:,7)/RF,depth,'--','Color','#4e5fa2','Parent',ax4,'linewidth',1.4);
cost = plot(radius_4_data(:,7)/RF,depth,'-','Color','#4e5fa2','Parent',ax4,'linewidth',1.4);
plot(base_case(7),base_case(1),'pentagram','MarkerFaceColor','#4e5fa2','MarkerEdgeColor','#4e5fa2','Parent',ax4);
set(ax4, 'YDir','reverse','Ylim', [min(depth),max(depth)],'xColor','#4e5fa2','gridcolor','k','xlim', [0,4500],'xtick',0:750:4500);
set(ax4,'YTickLabel',[]);
xlabel('SpCC [$ W_{e}^{-1}]', 'FontSize',11) % x-axis label
grid on

lgd = legend([con1 con2 con3 con4 star],"4-inch Well Diameter", "6-inch Well Diameter", "9^{5/8}-inch Well Diameter", "0.5 m Well Diameter", "Base Case",'numcolumns',3);
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 11;
lgd.FontSize = 11;
lgd.Position = [0.15 0.6 0.3 0.3]; %[110 500 150 90];
set(lgd,'color','white');

caption = {'{\bfProperties}'...
    ,'Massflow to Minimize Cost'...
    ,'Lateral Well Length to Minimize Cost'...
    ,'CO_{2} Working Fluid'...
    ,'35 �C km^{-1} Geothermal Gradient'...
    ,'15 �C Surface Temp'...
    ,'Four Horizontal Laterals'...
    ,strcat(['Year-' num2str(30) ' Depletion Values'])...
    ,'Malek|Adams|Rossi|Schiegg|Saar (2021)'};
annot = annotation('textbox','String',caption);
annot.Position = [0.05 0.68 0.3 0.3];
annot.FitBoxToText = 'on';
annot.BackgroundColor = 'White';
annot.FontSize = 11;

annota = annotation('textbox',[0.058 0.67 0.05 0.05],'String',"(a)",'EdgeColor','none','FontSize',14,'fontweight', 'bold');
annotb = annotation('textbox',[0.57 0.67 0.05 0.05],'String',"(b)",'EdgeColor','none','FontSize',14,'fontweight', 'bold');

% %% Specific Power Plot
% figure(2)
% four = plot(radius_1_data(:,5),depth,':','Color','#00B0F0','linewidth',1.4);
% hold on
% six = plot(radius_2_data(:,5),depth,'-.','Color','#00B0F0','linewidth',1.4);
% nine = plot(radius_3_data(:,5),depth,'--','Color','#00B0F0','linewidth',1.4);
% half = plot(radius_4_data(:,5),depth,'-','Color','#00B0F0','linewidth',1.4);
% base = plot(base_case(5),base_case(1),'pentagram','MarkerFaceColor','#00B0F0','MarkerEdgeColor','black');
% lgd = legend([four six nine half base], "4-inch Well Diameter", "6-inch Well Diameter", "9-5/8-inch Well Diameter", "0.5 m Well Diameter", "Base Case");
% lgd.Title.String = 'Legend';
% lgd.Title.FontSize = 11;
% lgd.FontSize = 11;
% lgd.Position = [0.15 0.6 0.3 0.3]; %[110 500 150 90];
% set(lgd,'color','white');
% set(gca, 'YDir','reverse','Ylim', [min(depth),max(depth)]);
% set(gca,'FontSize',11)
% xlabel('Specific Electric Power [W_{e}/m]', 'FontSize',11) % x-axis label
% ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
% title('Specific Power Variation for Different Radii', 'FontSize',11) % y-axis label
% grid on

saveas(figure(1),'images\radius.fig');
%saveas(figure(2),'images\radius_Specific_Power.fig');
