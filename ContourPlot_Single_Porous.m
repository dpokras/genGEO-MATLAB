clear all;
close all;

%fluid = 'Water';
fluid = 'CO2';
useBest = false;
Xtype = 'Transmissivity';
%Xtype = 'ResLength';
%Ztype = 'Power';
Ztype = 'LCOE';
%Ztype = 'SpCC';
dataPoints = 'Middleton';
%dataPoints = 'USGS';
%dataPoints = 'BaseCase';

%xlsxFile = 'data\Porous_Feb2020.xlsx';
xlsxFile = 'data\Porous_Jun2020.xlsx';
%xlsxFile = 'data\Conduction_Jun2020.xlsx';
co2Sheet = 'CO2';
waterSheet = 'Water (R245fa)';
%waterSheet = 'Water (R600a)';

[X,Y,Z,Z_Combined_Mask] = GetContourData(fluid, useBest, Xtype, Ztype, xlsxFile, co2Sheet, waterSheet);

if (strcmp(Ztype,'Power'))
    %Data is in MWe
    Z = Z*1000;
    contourLevels = [10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000];
elseif (strcmp(Ztype,'SpCC'))
    contourLevels = [1500, 2000, 3000, 3500, 4000, 5000, 6000, 8000, 10000, 15000, 20000, 50000];
elseif (strcmp(Ztype,'LCOE'))
    contourLevels = [50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 400, 500, 750, 1000, 2000];
end

% Contour Plot CO2
[M1,c1] = contour(X,Y,Z);
%[M1,c1] = contour(X,Y,Z_Combined);
%[M1,c1] = contour(X,Y,Z,30,'ShowText','on');
c1.LevelList = contourLevels;
c1.LineColor = 'k';
%c1.LineColor = '#F79646';
c1.LineWidth = 1;
%c1.LineStyle = '--';
c1.DisplayName = fluid;
if (strcmp(fluid,'CO2'))
    c1.DisplayName = 'CO_{2}';
end
clabel(M1,c1,contourLevels,'Margin',5,'LabelSpacing',400,'FontSize',12);


% Add dataPoints
hold on
fig = gcf;
legendLines = {};

%if (strcmp(dataPoints,'BaseCase'))
scatter(fig.CurrentAxes, [4.18], [2500], 100, 'pentagram', 'filled', 'm', 'DisplayName', 'Base-case');
legendLines = {'Base-case: (2.5km, 50mD, 300m)'};

if (strcmp(dataPoints,'USGS'))
    % Add USGS Data
    [usgs_data, usgs_headers, usgs_complete] = xlsread('data\USGS_Table_1.xlsx', 'Table_1');
    %thickness is column 4 (2 without row headers)
    thickness = usgs_data(:,2) * 0.3048;
    %perm is column 7 (really 5)
    perm = usgs_data(:,5);
    transmissivity = thickness .* perm;
    %depth is column 9 (really 7)
    depth = usgs_data(:,7) * 0.3048;
    scatter(fig.CurrentAxes, log10(transmissivity), depth, 'o', 'b', 'DisplayName', 'USGS (2012)');
    legendLines = {'DataPoints from USGS (2012)'};

elseif (strcmp(dataPoints,'Middleton'))
    %Add Other Data
    [x_data, x_headers, x_complete] = xlsread('data\transdata_all_data.xlsx', 'Sheet1');
    depth = x_data(:,1);
    transmissivity = x_data(:,2);
    scatter(fig.CurrentAxes, log10(transmissivity), depth, 'o', 'b', 'DisplayName', 'All-Middleton|LANL');
    [x_data, x_headers, x_complete] = xlsread('data\transdata_all_go.xlsx', 'Sheet1');
    depth = x_data(:,1);
    transmissivity = x_data(:,2);
    scatter(fig.CurrentAxes, log10(transmissivity), depth, 'o', 'r', 'DisplayName', 'Go-Middleton|LANL');
    legendLines = {'DataPoints from Middleton|LANL'};
end


FormatFigures(fluid, Xtype, Ztype, legendLines);


[~,dataFileName,~] = fileparts(xlsxFile);
fileName = strcat(['images\' Ztype '_Contour_' fluid '_' dataFileName '_' dataPoints '.png']);
saveas(gcf,fileName);

