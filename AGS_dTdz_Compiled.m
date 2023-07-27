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
dTdz_25_data = xlsread('data\dTdz_25.xlsx');
dTdz_30_data = xlsread('data\dTdz_30.xlsx');
dTdz_35_data = xlsread('data\dTdz_35.xlsx');
dTdz_40_data = xlsread('data\dTdz_40.xlsx');
% Data is in form: [Depth, optimum Res-length,optimum m-dot,Power, Specific Power, Prod_T, LCOE]
depth = dTdz_25_data(:,1)/1000;

%% Extract Data for Base Case
base_case = xlsread('data\Base Case Results.xlsx');
base_case(1) = base_case(1)/1000;
%% Plot Figures
figure(1)
subplot(2,2,3)
%reservoir legths
axis tight
plot(dTdz_25_data(:,2)/1000,depth,':k','linewidth',1.4);
hold on
plot(dTdz_30_data(:,2)/1000,depth,'-.k','linewidth',1.4);
plot(dTdz_35_data(:,2)/1000,depth,'--k','linewidth',1.4);
plot(dTdz_40_data(:,2)/1000,depth,'-k','linewidth',1.4);
plot(base_case(2)/1000,base_case(1),'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black');
xlabel('Optimal Lateral Well Length [km]', 'FontSize',10, 'HorizontalAlignment','right'); % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
set(gca,'FontSize',11,'xColor','k')
grid on
ax1 = gca;
set(ax1, 'Position', [0.04 0.09 0.40 0.56], 'YDir','reverse','Ylim', [min(depth),max(depth)],'Xlim',[0,15],'xtick',0:3:9,'gridcolor','k');
ax2 = axes('Position', [0.04 0.09 0.40 0.56],...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
set(ax2,'FontSize',11)
hold(ax2,'on')
%opt m-dot
dTdz_30_data(:,2)
plot(dTdz_25_data(:,3),depth,':','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(dTdz_30_data(:,3),depth,'-.','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(dTdz_35_data(:,3),depth,'--','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(dTdz_40_data(:,3),depth,'-','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(base_case(3),base_case(1),'pentagram','MarkerFaceColor','#4e5fa2','MarkerEdgeColor','#4e5fa2','Parent',ax2);
set(ax2, 'YDir','reverse','Ylim', [min(depth),max(depth)],'Xlim',[-100,150],'xColor','#4e5fa2','xtick',0:50:150,'gridcolor','k');
set(ax2,'YTickLabel',[]);
ylim([1,8])
xlh = xlabel('Optimal Mass Flowrate [kg s^{-1}]', 'FontSize',11, 'HorizontalAlignment','left'); % x-axis label
grid on

% Power Subplot
subplot(2,2,[2,4])
con1 = plot(NaN,NaN,':k','linewidth',1.4);
hold on
con2 = plot(NaN,NaN,'-.k','linewidth',1.4);
con3 = plot(NaN,NaN,'--k','linewidth',1.4);
con4 = plot(NaN,NaN,'-k','linewidth',1.4);
space = plot(NaN,NaN,'-w','linewidth',1.4);
star = plot(NaN,NaN,'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black');
optreslength = plot(NaN,NaN,'-m','linewidth',1.4);
optmdot = plot(NaN,NaN,'-b','linewidth',1.4);
plot(dTdz_25_data(:,4)/2,depth,':','Color','black','linewidth',1.4);
plot(dTdz_30_data(:,4)/2,depth,'-.','Color','black','linewidth',1.4);
plot(dTdz_35_data(:,4)/2,depth,'--','Color','black','linewidth',1.4);
pow = plot(dTdz_40_data(:,4)/2,depth,'-','Color','black','linewidth',1.4);
plot(base_case(4),base_case(1),'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black')
set(gca,  'Position', [0.57 0.09 0.40 0.56], 'YDir','reverse','Ylim', [min(depth),max(depth)],'xlim',[0,4],'xColor','black','xtick',0:0.8:4,'gridcolor','k');
grid on
set(gca,'FontSize',11)
xlabel('Electric Power [MW_{e}]', 'FontSize',11) % x-axis label
ylabel('Vertical Well Depth [km]', 'FontSize',11) % y-axis label
ax1 = gca;
ax1_pos = ax1.Position;
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
set(ax2,'FontSize',11)
hold(ax2,'on')
% SpCC
plot(dTdz_25_data(:,7)/RF,depth,':','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(dTdz_30_data(:,7)/RF,depth,'-.','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(dTdz_35_data(:,7)/RF,depth,'--','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
cost = plot(dTdz_40_data(:,7)/RF,depth,'-','Color','#4e5fa2','Parent',ax2,'linewidth',1.4);
plot(base_case(7),base_case(1),'pentagram','MarkerFaceColor','#4e5fa2','MarkerEdgeColor','#4e5fa2','Parent',ax2)
set(ax2, 'YDir','reverse','Ylim', [min(depth),max(depth)],'xColor','#4e5fa2','gridcolor','k');
set(ax2,'YTickLabel',[]);
xlabel('SpCC [$ W_{e}^{-1}]', 'FontSize',11) % x-axis label
xlim([0,5000])
grid on

lgd = legend([con1 con2 con3 con4 star],"\nabla T = 25 \circC km^{-1}", "\nabla T = 30 \circC km^{-1}", "\nabla T = 35 \circC km^{-1}", "\nabla T = 40 \circC km^{-1}", "Base Case",'numcolumns',3);
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 11;
lgd.FontSize = 11;
lgd.Position = [0.15 0.6 0.3 0.3]; %[110 500 150 90];
set(lgd,'color','white');

caption = {'{\bfProperties}'...
    ,'Massflow to Minimize Cost'...
    ,'Lateral Well Length to Minimize Cost'...
    ,'CO_{2} Working Fluid'...
    ,'15 �C Surface Temp'...
    ,'0.5 m Well Diameter'...
    ,'Four Horizontal Laterals'...
    ,strcat(['Year-' num2str(30) ' Depletion Values'])...
    ,'Malek|Adams|Rossi|Schiegg|Saar (2021)'};
annot = annotation('textbox','String',caption);
annot.Position = [0.05 0.68 0.3 0.3];
annot.FitBoxToText = 'on';
annot.BackgroundColor = 'White';
annot.FontSize = 11;

annota = annotation('textbox',[0.045 0.67 0.05 0.05],'String',"(a)",'EdgeColor','none','FontSize',14,'fontweight', 'bold');
annotb = annotation('textbox',[0.57 0.67 0.05 0.05],'String',"(b)",'EdgeColor','none','FontSize',14,'fontweight', 'bold');

saveas(figure(1),'images\dTdz.fig');