function [WellCost] = CapitalCost_Well(welldepth, wellLength, wellradius, params)

    PPI_current = PPIs.return_PPI('PPI_O&G', params.costYear, params); %2022 2.195
    
    if welldepth > 0
        if (strcmp(params.wellCostType,'Ideal'))
            % revised GETEM cost correllations 2022 large diameter -
            % 12.25in vertical well)
            PPI_ref = PPIs.return_PPI('PPI_O&G', 2003, params); %2003 153.492 
            C_well = PPI_current/PPI_ref * (-62.2.*welldepth + 1290.*(2.*wellradius).*welldepth + 275300);
            %disp('Caution: Using Ideal Well Relations');
        elseif (strcmp(params.wellCostType,'Baseline'))
            PPI_ref = PPIs.return_PPI('PPI_O&G', 2003, params); %2003 153.492 
%             C_well =   PPI_O_G .* (0.40.*(welldepth).^2 + 203.95.*welldepth + 1,589,396.86)
            C_well = PPI_current/PPI_ref .* (0.105.*(welldepth).^2 + 1776.*(2.*wellradius).*welldepth + 275300);
        else
            throw(MException('CapitalCost_Well:unknownWellType','Unknown Well Type'));
        end
        dC_well_CO2 = PPI_current/PPI_ref .* (265.*(2.*wellradius)*welldepth + 133.*welldepth);
    else
        % horizontal well (revised GETEM correllation 2022 small diameter -
        % 8.5 in horizontal well)

        PPI_ref = PPIs.return_PPI('PPI_O&G', 2003, params); %2003 153.492 
        %horizontal well cost modifier = 15%
        X_HWell = 1.15;
        C_well = PPI_current./PPI_ref .* X_HWell.* (0.105.*(wellLength).^2 + 1776.*(2.*wellradius).*wellLength + 275300);
%         PPI_ref = PPI('PPI_O&G', 2022); %2003 153.492 
%         C_well =   PPI_current/PPI_ref  .* 0.10.*(wellLength).^2 + 1280.*wellLength + 273500;
        dC_well_CO2 =   PPI_current/PPI_ref .* (265.*(2.*wellradius).*wellLength + 133.*wellLength);
    end
    if (strcmp(params.fluid,'Water'))
        WellCost = C_well ./ params.SuccessRate_well;
    elseif (strcmp(params.fluid,'CO2'))
        WellCost = (C_well + dC_well_CO2) ./ params.SuccessRate_well;
    else
        throw(MException('CapitalCost_Wellfield:NotImplemented','Not Implemented'));
    end
end