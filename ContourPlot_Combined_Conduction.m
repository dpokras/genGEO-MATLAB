clear all;
close all;

%Xtype = 'Transmissivity';
Xtype = 'ResLength';

%xlsxFile = 'data\Conduction_Jun2020.xlsx';
xlsxFile = 'data\Conduction_Oct2020.xlsx';
co2Sheet = 'CO2 (Baseline)';
waterSheet = 'Water (Baseline)';
%co2Sheet = 'CO2 (Ideal)';
%waterSheet = 'Water (Ideal)';

% Define plot limits
XLim = 10000;
YLim = 8000;
WACC = 9.6;
F_OM = 4.5;
year = 30;

for i=1:2
    
    if (i==1)
        Ztype = 'Power';
    elseif (i==2)
        Ztype = 'LCOE_Greenfield';
    elseif (i==3)
        Ztype = 'SpCC_Greenfield';
    end

    figure(i);
    
    [X_Water,Y_Water,Z_Water,Z_Combined_Mask_Water] = GetContourData('Water', false, Xtype, Ztype, xlsxFile, co2Sheet, waterSheet);
    [X_CO2,Y_CO2,Z_CO2,Z_Combined_Mask_CO2] = GetContourData('CO2', false, Xtype, Ztype, xlsxFile, co2Sheet, waterSheet);

    if (strcmp(Ztype,'Power'))
        %Data is in MWe
        Z_CO2 = Z_CO2*1000;
        Z_Water = Z_Water*1000;
        contourLevels = [10, 32, 100, 320, 1000, 3200, 10000, 18000, 32000, 56000, 100000];
    elseif (contains(Ztype,'SpCC'))
        contourLevels = [1500, 2000, 3000, 3500, 4000, 5000, 6000, 8000, 10000, 15000, 20000, 35000, 50000, 75000, 1e5, 1.5e5 2e5, 3e5];
    elseif (contains(Ztype,'LCOE'))
        contourLevels = [50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 4000, 7000, 10000, 20000];
    end

    % Contour Plot CO2
    [M1,c1] = contour(X_CO2,Y_CO2,Z_CO2);
    c1.LevelList = contourLevels;
    c1.LineColor = '#F79646';
    c1.LineWidth = 1;
    %c1.LineStyle = '--';
    c1.DisplayName = 'CO_{2}';
    clabel(M1,c1,contourLevels,'Margin',5,'LabelSpacing',400,'FontSize',12);

    % Contour Plot Water
    hold on
    [M2,c2] = contour(X_Water,Y_Water,Z_Water);
    c2.LevelList = contourLevels;
    c2.LineColor = 'k';
    c2.LineWidth = 1;
    c2.DisplayName = 'Water';
    clabel(M2,c2,contourLevels,'Margin',5,'LabelSpacing',400,'FontSize',12);

    if (contains(co2Sheet,'Ideal'))
        legendLines = {'Ideal Well Cost'};
    else
        legendLines = {'Baseline Well Cost'};
    end

    FormatFigures('CO2&Water', Xtype, Ztype, legendLines, XLim, YLim, year, WACC, F_OM);


    [~,dataFileName,~] = fileparts(xlsxFile);
    if (contains(co2Sheet,'Ideal'))
        fileName = strcat(['images\' Ztype '_Contour_Ideal_' dataFileName '.png']);
    else
        fileName = strcat(['images\' Ztype '_Contour_Baseline_' dataFileName '.png']);
    end
    saveas(gcf,fileName);

end
