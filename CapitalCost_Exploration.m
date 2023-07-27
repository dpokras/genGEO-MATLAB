function [C_exploration] = CapitalCost_Exploration(params)

    PPI_OG_s = PPI('PPI_O&G-s', params.costYear); %2019 1.437
    X_PC_expl = 1.15;
    X_IC_expl = 1.05;

%     %WellFieldLength is the length of the wellfield, assuming it's a square
%     if (strcmp(params.wellFieldType,'Doublet'))
%         WellFieldLength_km = 1;
%     elseif (strcmp(params.wellFieldType,'5spot'))
%         WellFieldLength_km = 1;
%     elseif (strcmp(params.wellFieldType,'5spot_SharedNeighbor'))
%         WellFieldLength_km = 1;
%     elseif (strcmp(params.wellFieldType,'5spot_ManyN'))
%         WellFieldLength_km = N_5spot;
%     elseif (strcmp(params.wellFieldType,'Tungsten'))
%         WellFieldLength_km = 1;
%     else
%         throw(MException('CapitalCost_Wellfield:BadWellFieldType','Bad WellFieldType'));
%     end

    
    G_characterizationWells = 2;

    L_wellZone = WellFieldLength_km * 1000;
    A_CO2_AMA = (L_wellZone + 1600)^2;

    C_modeling = X_IC_expl * X_PC_expl * PPI_OG_s * 508000;
    C_well = CapitalCost_Well(params.depth, 2*params.well_radius, 'Water', 1, params);
    dC_characterizationWells = 0.2 * C_well * G_characterizationWells / params.SuccessRate_well;
    dC_modeling_CO2 = X_IC_expl * X_PC_expl * PPI_OG_s * 44800 * (1/1e6) * A_CO2_AMA;


    if (strcmp(params.fluid,'Water'))
        C_exploration = C_modeling + dC_characterizationWells;
    elseif (strcmp(params.fluid,'CO2'))
        C_exploration = C_modeling + dC_characterizationWells + dC_modeling_CO2;
    else
        throw(MException('CapitalCost_Exploration:NotImplemented','Not Implemented'));
    end

end