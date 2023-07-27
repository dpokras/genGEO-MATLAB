function result = CapitalCost_Water(Q_preheater_total, Q_boiler_total, W_turbine_total, Q_recuperator_total, Q_desuperheater_total, ...
    Q_condenser_total, W_pump_orc_total, W_pump_prod_total, ...
    dT_range_CT, dT_LMTD_preheater, dT_LMTD_boiler, dT_LMTD_recuperator, params)

    % All heat/power in Watts

    %Check that simulation is for water
    if (strcmp(params.fluid,'Water') == false)
        throw(MException('CapitalCost_Water:Wrong_Fluid','Wrong Fluid--Not Water'));
    end


    % Capital Cost is
    % 1. Surface Plant
    % 2. Gathering System
    % 3. Wells
    % 4. Wellfield Development
    % 5. Exploration
    % 6. Well Stimulation


    % 1. Surface Plant
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    result_surfacePlant = CapitalCost_SurfacePlant_ORC(Q_preheater_total, Q_boiler_total, W_turbine_total, Q_recuperator_total, Q_desuperheater_total, Q_condenser_total, ...
        W_pump_orc_total, W_pump_prod_total, dT_range_CT, dT_LMTD_preheater, dT_LMTD_boiler, dT_LMTD_recuperator, params);

    % 2. Gathering System
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (params.hasSurfaceGatheringSystem == true)
        C_gatheringSystem = CapitalCost_SurfacePipe(params);
    else
        C_gatheringSystem = 0;
    end

    % 3. Well Cost
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Number of wells to pay for
    if (strcmp(params.wellFieldType,'Doublet'))
        G_wells_injection = 1;
        G_wells_production = 1;
    elseif (strcmp(params.wellFieldType,'5spot'))
        G_wells_injection = 1;
        G_wells_production = 4;
    elseif (strcmp(params.wellFieldType,'5spot_SharedNeighbor'))
        G_wells_injection = 1;
        G_wells_production = 1;
    elseif (strcmp(params.wellFieldType,'5spot_ManyN'))
        G_wells_injection = N_5spot^2;
        G_wells_production = (N_5spot+1)^2;
    elseif (strcmp(params.wellFieldType,'Tungsten'))
        G_wells_injection = 4;
        G_wells_production = 4;
    else
        throw(MException('CapitalCost_Water:BadWellFieldType','Bad WellFieldType'));
    end

    %If conduction system, the wellLength is the depth plus half reservoir
    %length
    C_well_vertical = CapitalCost_Well(params.depth, 2*params.well_radius, params.fluid, params.SuccessRate_well, params);
    C_well_horizontal = CapitalCost_Well((params.depth+0.5*params.res_length), 2*params.well_radius, params.fluid, params.SuccessRate_well, params) - C_well_vertical;
    
    if (strcmp(params.system,'Porous'))
        C_well = C_well_vertical;
    else
        C_well = C_well_vertical + params.n_streams*C_well_horizontal
%     elseif (strcmp(params.system,'Conduction1'))
%         C_well = C_well_vertical + C_well_horizontal;
%     elseif (strcmp(params.system,'Conduction2'))
%         C_well = C_well_vertical + 2*C_well_horizontal;
%     elseif (strcmp(params.system,'Conduction4'))
%         C_well = C_well_vertical + 4*C_well_horizontal;
%     elseif (strcmp(params.system,'Conduction8'))
%         C_well = C_well_vertical + 8*C_well_horizontal;
%     else
%         throw(MException('CapitalCost_Water:NotImplemented','CapitalCost_Water:Not Implemented'));
    end

    C_wells_injection =  G_wells_injection * C_well;
    C_wells_production = G_wells_production * C_well;

    % 4. Wellfield Development 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C_wellfield = CapitalCost_Wellfield(params);

    % 5. Exploration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % No exploration costs if conduction system
    if (strcmp(params.system,'Porous'))
        C_exploration = CapitalCost_Exploration(params);
    elseif (contains(params.system,'Conduction'))
        C_exploration = 0;
    else
        throw(MException('CapitalCost_Water:NotImplemented','Not Implemented'));
    end

    % 6. Well Stimulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C_stimulation = 0;


    % SUM IT UP
    C_brownfield = result_surfacePlant.C_plant + C_gatheringSystem + C_wells_production;
    C_greenfield = result_surfacePlant.C_plant + C_gatheringSystem + C_wells_production + C_wells_injection + C_wellfield + C_exploration + C_stimulation;

    result.CostSurfacePlant = result_surfacePlant;
    result.C_surfacePlant = result_surfacePlant.C_plant;
    result.C_gatheringSystem = C_gatheringSystem;
    result.C_wells_production = C_wells_production;
    result.C_wells_injection = C_wells_injection;
    result.C_wellfield = C_wellfield;
    result.C_exploration = C_exploration;
    result.C_stimulation = C_stimulation;
    result.C_brownfield = C_brownfield;
    result.C_greenfield = C_greenfield;

end

