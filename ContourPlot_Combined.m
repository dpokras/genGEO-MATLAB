clear all;
close all;

Xtype = 'Transmissivity';
%Xtype = 'ResLength';
%Ztype = 'Power';
Ztype = 'LCOE';
%Ztype = 'SpCC_Brownfield';
%Ztype = 'SpCC_Greenfield';
%dataPoints = 'Middleton';
dataPoints = 'USGS';
%dataPoints = 'BaseCase';
%dataPoints = 'NorthSea';

%xlsxFile = 'data\Porous_Feb2020.xlsx';
%xlsxFile = 'data\Porous_Jun2020.xlsx';
%xlsxFile = 'data\Porous_Jul2020.xlsx';
%xlsxFile = 'data\Porous_Aug2020.xlsx';
%xlsxFile = 'data\Porous_dTdz30_Oct2020.xlsx';
xlsxFile = 'data\Porous_dTdz35_Oct2020.xlsx';
%xlsxFile = 'data\Porous_dTdz40_Oct2020.xlsx';
%xlsxFile = 'data\Conduction_Jun2020.xlsx';
co2Sheet = 'CO2';
waterSheet = 'Water (R245fa)';
%waterSheet = 'Water (R600a)';

% Define plot limits
%XLim = 5;
XLim = 7;
YLim = 7000;
%WACC = 9.6;
WACC = 4;
%F_OM = 4.5;
F_OM = 5.5;
dTdz = 35;
year = 1;
showProperties = true;

% LCOE_factor to adjust to 4% and 7%
if (WACC == 9.6 && F_OM == 4.5)
    LCOE_factor = 1;
elseif (WACC == 4 && F_OM == 5.5)
    LCOE_factor = 0.701;
else
    throw(MException('ContourPlot_Combined:NotImplemented','Not Implemented'));
end

if (contains(Ztype,'LCOE'))
    % Find CO2 Brownfield/Greenfield and only water greenfield
    [X_Water,Y_Water,Z_Water,Z_Combined_Mask_Water] = GetContourData('Water', false, Xtype, 'LCOE_Greenfield', xlsxFile, co2Sheet, waterSheet);
    [X_CO2_G,Y_CO2_G,Z_CO2_G,Z_Combined_Mask_CO2_G] = GetContourData('CO2', false, Xtype, 'LCOE_Greenfield', xlsxFile, co2Sheet, waterSheet);
    [X_CO2,Y_CO2,Z_CO2,Z_Combined_Mask_CO2] = GetContourData('CO2', false, Xtype, 'LCOE_Brownfield', xlsxFile, co2Sheet, waterSheet);
else
    [X_Water,Y_Water,Z_Water,Z_Combined_Mask_Water] = GetContourData('Water', false, Xtype, Ztype, xlsxFile, co2Sheet, waterSheet);
    [X_CO2,Y_CO2,Z_CO2,Z_Combined_Mask_CO2] = GetContourData('CO2', false, Xtype, Ztype, xlsxFile, co2Sheet, waterSheet);
end


if (strcmp(Ztype,'Power'))
    %Data is in MWe
    Z_CO2 = Z_CO2*1000;
    Z_Water = Z_Water*1000;
    contourLevels = [10, 32, 100, 320, 1000, 3200, 10000, 18000, 32000, 56000, 100000];
elseif (contains(Ztype,'SpCC'))
    contourLevels = [1500, 2000, 3000, 3500, 4000, 5000, 6000, 8000, 10000, 15000, 20000, 50000];
elseif (contains(Ztype,'LCOE'))
    Z_CO2 = LCOE_factor * Z_CO2;
    Z_CO2_G = LCOE_factor * Z_CO2_G;
    Z_Water = LCOE_factor * Z_Water;
    contourLevels = [35, 40, 50, 75, 100, 125, 200, 500, 2000];
end

% Contour Plot CO2

[M1,c1] = contour(X_CO2,Y_CO2,Z_CO2);
c1.LevelList = contourLevels;
c1.LineColor = '#F79646';
c1.LineWidth = 1;
c1.DisplayName = 'CO_2';
clabel(M1,c1,contourLevels,'Margin',5,'LabelSpacing',400,'FontSize',12);
hold on;

% Add CO2 greenfield, if LCOE selected
if (contains(Ztype,'LCOE'))
    [M2,c2] = contour(X_CO2_G,Y_CO2_G,Z_CO2_G);
    c2.LevelList = [75, 100, 200, 500];
    c2.LineColor = '#F79646';
    c2.LineWidth = 0.5;
    c2.LineStyle = '-.';
    c2.DisplayName = 'CO_2-GF';
    clabel(M2,c2,contourLevels,'Margin',5,'LabelSpacing',400,'FontSize',12);
end

% Contour Plot Water
[M3,c3] = contour(X_Water,Y_Water,Z_Water);
c3.LevelList = contourLevels;
c3.LineColor = 'k';
c3.LineWidth = 1;
c3.DisplayName = 'Water';
clabel(M3,c3,contourLevels,'Margin',5,'LabelSpacing',400,'FontSize',12);

% Add Threshold
if (size(Z_Water,1)==size(Z_CO2,1) && size(Z_Water,2)==size(Z_CO2,2))
    % Set values less than or equal to zero to NaN
    Z_Water(Z_Water <= 0) = NaN;
    Z_CO2(Z_CO2 <= 0) = NaN;
    Z_diff = Z_Water - Z_CO2;
    [M4,c4] = contour(X_Water,Y_Water,Z_diff);
    c4.LevelList = [0];
    c4.LineColor = 'r';
    c4.LineWidth = 2.4;
    c4.LineStyle = ':';
    c4.DisplayName = 'Threshold';
    % if LCOE, also plot GF threshold
    if (contains(Ztype,'LCOE'))
        Z_CO2_G(Z_CO2_G <= 0) = NaN;
        Z_diff_G = Z_Water - Z_CO2_G;
        [M5,c5] = contour(X_Water,Y_Water,Z_diff_G);
        c5.LevelList = [0];
        c5.LineColor = 'r';
        c5.LineWidth = 2;
        c5.LineStyle = '-.';
        c5.DisplayName = 'GF Threshold';
    end
end

% Add datapoints to plot
fig = gcf;
legendLines = {};

%if (strcmp(dataPoints,'BaseCase'))
scatter(fig.CurrentAxes, [4.18], [2500], 100, 'pentagram', 'filled', 'm', 'DisplayName', 'Base-case');
if (contains(waterSheet,'R245fa'))
    legendLines = {'Base-case: (2.5km, 15000mD-m)', 'Water System ORC Fluid: R245fa'};
elseif (contains(waterSheet,'R600a'))
    legendLines = {'Base-case: (2.5km, 15000mD-m)', 'Water System ORC Fluid: R600a'};
end

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
    z_q = interp2(X_CO2,Y_CO2,Z_CO2,log10(transmissivity),depth);
    scatter(fig.CurrentAxes, log10(transmissivity), depth, 'o', 'b', 'DisplayName', 'USGS (2013)');

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
    
elseif (strcmp(dataPoints,'NorthSea'))
    %Add Other Data
    [x_data, x_headers, x_complete] = xlsread('data\N_Sea_CO2_geothermal_data_gluyas_v2.xlsx', 'Sheet1');
    % EISB area is rows 1 to 6
    depth_EISB = x_data(1:6,1);
    thickness_EISB = x_data(1:6,2);
    perm_EISB = x_data(1:6,3);
    transmissivity_EISB = thickness_EISB .* perm_EISB;
    scatter(fig.CurrentAxes, log10(transmissivity_EISB), depth_EISB, 'o', 'MarkerEdgeColor', '#0072BD', 'DisplayName', 'EISB-NorthSea|Gluyas');
    % Central area is rows 7 to 56
    depth_Central = x_data(7:56,1);
    thickness_Central = x_data(7:56,2);
    perm_Central = x_data(7:56,3);
    transmissivity_Central = thickness_Central .* perm_Central;
    scatter(fig.CurrentAxes, log10(transmissivity_Central), depth_Central, 'o', 'MarkerEdgeColor', '#A2142F', 'DisplayName', 'Central-NorthSea|Gluyas');
    % EISB area is rows 57 to 69
    depth_Southern = x_data(57:69,1);
    thickness_Southern = x_data(57:69,2);
    perm_Southern = x_data(57:69,3);
    transmissivity_Southern = thickness_Southern .* perm_Southern;
    scatter(fig.CurrentAxes, log10(transmissivity_Southern), depth_Southern, 'o', 'MarkerEdgeColor', '#77AC30', 'DisplayName', 'Southern-NorthSea|Gluyas');
end

FormatFigures('CO2&Water', Xtype, Ztype, legendLines, XLim, YLim, year, WACC, F_OM, dTdz, showProperties);


[~,dataFileName,~] = fileparts(xlsxFile);
fileName = strcat(['images\' Ztype '_Contour_Overlap_' dataFileName '_' dataPoints '.png']);
saveas(gcf,fileName);

