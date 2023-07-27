function [caption, title] = FigureCaption(Ztype, fluid, legendLines, year, WACC, F_OM, dTdz)

if (contains(Ztype,'Brownfield'))
    prefix = 'Brownfield ';
elseif (contains(Ztype,'Greenfield'))
    prefix = 'Greenfield ';
else
    prefix = '';
end

if (WACC == 9.6 && F_OM == 4.5)
    lcoeType = 'LCOE_{LAZARD}';
elseif (WACC == 4 && F_OM == 5.5)
    lcoeType = 'LCOE_{LIKELY}';
else
    lcoeType = 'LCOE';
end


if (strcmp(Ztype,'Power'))
    if (strcmp(fluid,'CO2'))
        title = 'Electricity Generated per 5-spot for CO_2 Geothermal Systems [kW\fontsize{12}e\fontsize{18}]';
    elseif (strcmp(fluid,'Water'))
        title = 'Electricity Generated per 5-spot for Water Geothermal Systems [kW\fontsize{12}e\fontsize{18}]';
    elseif (strcmp(fluid,'CO2&Water'))
        title = 'Electricity Generated per 5-spot for CO_2 and Water Geothermal Systems [kW\fontsize{12}e\fontsize{18}]';
    end
    caption = {'{\bfProperties}'...
        ,strcat([ num2str(dTdz) ' °C/km Gradient, 15°C Surface']) ...
        ,'41 cm Well Diameter'...
        ,'Massflow to Maximize Power'...
        ,'5-spot Well Spacing (1 km²)'...
        ,strcat(['Year-' num2str(year) ' Depletion Values'])...
        ,'Adams|Saar, ETH-Zürich (geg.ethz.ch)'};
elseif (contains(Ztype,'SpCC'))
    if (strcmp(fluid,'CO2'))
        title = strcat([ prefix 'Specific Capital Cost for CO_2 Geothermal Systems [2019$/kW\fontsize{12}e\fontsize{18}]']);
    elseif (strcmp(fluid,'Water'))
        title = strcat([ prefix 'Specific Capital Cost for Water Geothermal Systems [2019$/kW\fontsize{12}e\fontsize{18}]']);
    elseif (strcmp(fluid,'CO2&Water'))
        title = strcat([ prefix 'Specific Capital Cost for CO_2 and Water Geothermal Systems [2019$/kW\fontsize{12}e\fontsize{18}]']);
    end
    caption = {'{\bfProperties}'...
  ,strcat([ num2str(dTdz) ' °C/km Gradient, 15°C Surface']) ...
        ,strcat([ num2str(dTdz) ' °C/km Gradient, 15°C Surface']) ...
        ,'41 cm Well Diameter'...
        ,'5-spot\_SN Well Spacing & Cost (1 km²)'...
        ,'Massflow to Minimize Cost'...
        ,strcat(['Year-' num2str(year) ' Depletion Values'])...
        ,'Adams|Saar, ETH-Zürich (geg.ethz.ch)'};
elseif (contains(Ztype,'LCOE'))
    if (strcmp(fluid,'CO2'))
        title = strcat([ prefix lcoeType ' for CO_2 Geothermal Systems [2019$/(MW\fontsize{12}e\fontsize{18}-h)]']);
    elseif (strcmp(fluid,'Water'))
        title = strcat([ prefix lcoeType ' for Water Geothermal Systems [2019$/(MW\fontsize{12}e\fontsize{18}-h)]']);
    elseif (strcmp(fluid,'CO2&Water'))
        title = strcat([ prefix lcoeType ' for CO_2 and Water Geothermal Systems [2019$/(MW\fontsize{12}e\fontsize{18}-h)]']);
    end
    caption = {'{\bfProperties}'...
        ,strcat([ num2str(dTdz) ' °C/km Gradient, 15°C Surface']) ...
        ,'41 cm Well Diameter'...
        ,'5-spot\_SN Well Spacing & Cost (1 km²)'...
        ,'Massflow to Minimize Cost'...
        ,strcat([ num2str(WACC) '% WACC, ' num2str(F_OM) '% O&M, 25-year Lifetime'])...
        ,strcat(['Year-' num2str(year) ' Depletion Values'])...
        ,'Adams|Saar, ETH-Zürich (geg.ethz.ch)'};
end



if (~isempty(legendLines))
    caption = vertcat(caption{1:4},legendLines(:),caption{5:end});
end
end

