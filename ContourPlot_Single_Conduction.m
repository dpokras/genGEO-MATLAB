clear all;
close all;

%fluid = 'Water';
fluid = 'CO2';
useBest = false;
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

for i = 1:2
    if (i == 1)
        Ztype = 'Power';
    elseif (i == 2)
        Ztype = 'LCOE_Greenfield';
    elseif (i == 3)
        Ztype = 'SpCC_Greenfield';
    end
    
    figure(i);
        
    [X,Y,Z,Z_Combined_Mask] = GetContourData(fluid, useBest, Xtype, Ztype, xlsxFile, co2Sheet, waterSheet);

    if (strcmp(Ztype,'Power'))
        %Data is in MWe
        Z = Z*1000;
        contourLevels = [10, 18, 32, 56, 100, 180, 320, 560, 1000, 1800, 3200, 5600, 10000, 18000, 32000, 56000, 100000];
    elseif (contains(Ztype,'SpCC'))
        contourLevels = [1500, 2000, 3000, 3500, 4000, 5000, 6000, 8000, 10000, 15000, 20000, 35000, 50000, 75000, 1e5, 1.5e5 2e5, 3e5];
    elseif (contains(Ztype,'LCOE'))
        contourLevels = [50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000, 8000, 10000];
    end

    % Contour Plot CO2
    [M1,c1] = contour(X,Y,Z);
    %[M1,c1] = contour(X,Y,Z_Combined);
    %[M1,c1] = contour(X,Y,Z,30,'ShowText','on');
    c1.LevelList = contourLevels;
    c1.LineColor = 'k';
    c1.LineWidth = 1;
    %c1.LineStyle = '--';
    clabel(M1,c1,contourLevels,'Margin',5,'LabelSpacing',400,'FontSize',12);

    if (strcmp(fluid,'CO2'))
        sheetName = regexprep(co2Sheet, ' ', '');
        if (contains(co2Sheet,'Ideal'))
            legendLines = {'Ideal Well Cost'};
        else
            legendLines = {'Baseline Well Cost'};
        end
    else
        sheetName = regexprep(waterSheet, ' ', '');
        if (contains(waterSheet,'Ideal'))
            legendLines = {'Ideal Well Cost'};
        else
            legendLines = {'Baseline Well Cost'};
        end
    end
    
    FormatFigures(fluid, Xtype, Ztype, legendLines, XLim, YLim, year, WACC, F_OM)
    
    
    [~,dataFileName,~] = fileparts(xlsxFile);
    fileName = strcat(['images/' Ztype '_Contour_' sheetName '_' dataFileName '.png']);
    saveas(gcf,fileName);
end

