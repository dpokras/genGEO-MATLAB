function result = CapitalCost_CPG(result_turbine, result_HEX, result_pumps, result_cooling_tower, result_tank, result_comp_unit, params)
%Power values are all in kW for a single Injection/Production (IP) pair.

    % Capital Cost is
    % 1. Surface Plant
    % 2. Gathering System
    % 3. Wells
    % 4. Wellfield Development
    % 5. Exploration
    % 6. Well Stimulation


    % 1. Surface Plant
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C_TBM_plant = CapitalCost_SurfacePlant_CPG(result_turbine, result_HEX, result_pumps, result_cooling_tower, result_tank, result_comp_unit, params);
    

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

    %If conduction system, the wellLength is the depth plus half reservoir
    %length
    C_well_vertical = CapitalCost_Well(params.depth, 0, mean([params.well_radius,params.side_stream_radius]), params);
    C_well_horizontal = CapitalCost_Well(0, params.res_length/2, params.side_stream_radius, params);

    C_well_vertical_green_saved = CapitalCost_Well(3000, 0, params.well_radius, params);
    C_well_vertical_grean_extended = CapitalCost_Well(params.depth(end), 0, params.well_radius, params) - C_well_vertical_green_saved;

    C_wells_horizontal = 2 * C_well_horizontal * params.n_streams;
    
    C_well_brownfield = C_well_vertical + C_well_vertical_grean_extended + C_wells_horizontal;
    C_well_greenfield = 2*C_well_vertical + C_wells_horizontal;

    % 4. Wellfield Development 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C_wellfield = CapitalCost_Wellfield(params);
    
    % 5. Exploration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % No exploration costs if conduction system

    C_exploration = 0;

    % 6. Well Stimulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C_stimulation = 0;

    % 7. CO2 Purcurement Costs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Liquid CO2 procurement costs
    C_tCO2 = 390 * result_tank.V_total_surface + result_tank.V_store; %2022 USD/MT 

    C_TBM_subsurface_brownfield = C_well_brownfield + C_exploration + C_stimulation;
    C_TBM_subsurface_greenfield = C_well_greenfield + C_wellfield + C_exploration + C_stimulation;
    C_TBM_brownfield = C_TBM_plant + C_gatheringSystem + C_TBM_subsurface_brownfield + C_tCO2;
    C_TBM_greenfield = C_TBM_plant + C_gatheringSystem + C_TBM_subsurface_greenfield + C_tCO2;
        
    % Calculate Total Captial Investment (TCI)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C_TCI = CapitalCost_SurfacePlantfunc.TCI([C_TBM_brownfield, C_TBM_greenfield]);
    C_TCI_SurfacePlant = CapitalCost_SurfacePlantfunc.TCI(C_TBM_plant+C_tCO2);
    C_TCI_Subsurface = CapitalCost_SurfacePlantfunc.TCI([C_TBM_subsurface_brownfield, C_TBM_subsurface_greenfield]);

    result.CostSurfacePlant = C_TCI_SurfacePlant;
    result.CostSubsurface_brownfield = C_TCI_Subsurface(1);
    result.CostSubsurface_greenfield = C_TCI_Subsurface(2);
    result.C_gatheringSystem = C_gatheringSystem;
    result.C_wells_production = C_well_vertical;
    result.C_wells_horizontal = C_wells_horizontal;
    result.C_wellfield = C_wellfield;
    result.C_exploration = C_exploration;
    result.C_stimulation = C_stimulation;
    result.C_brownfield = C_TCI(1);
    result.C_greenfield = C_TCI(2);

end

