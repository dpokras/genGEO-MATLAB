function result = total_analytic_system_GEOPHIRES(params)

    %Check that simulation is for water
    if (strcmp(params.fluid,'Water') == false)
        throw(MException('total_analytic_system_GEOPHIRES:Wrong_Fluid','Wrong Fluid--Not Water'));
    end

        %Prod/Inj Well Multipliers
    if (strcmp(params.wellFieldType, 'Doublet'))
        Field_IP_Multiplier = 1;
        numInjWells = 1;
        numProdWells = 1;
    elseif (strcmp(params.wellFieldType, '5spot_SharedNeighbor'))
        Field_IP_Multiplier = 4;
        numInjWells = 1;
        numProdWells = 1;
    elseif (strcmp(params.wellFieldType, '5spot'))
        Field_IP_Multiplier = 4;
        numInjWells = 1;
        numProdWells = 1;
    elseif (strcmp(params.wellFieldType, '5spot_ManyN'))
        Field_IP_Multiplier = 4 * N_5spot^2;
        numInjWells = 1;
        numProdWells = 1;
    elseif (strcmp(params.wellFieldType, 'Tungsten'))
        Field_IP_Multiplier = 1;
        numInjWells = 4;
        numProdWells = 4;
    else
        throw(MException('total_analytic_system_water:BadWellfieldType','Bad Wellfield Type'));
    end
    
    % Get reservoir Impedance
    P_f_initial = params.depth * 1000 * 9.81;
    %T_f_initial = params.T_surface_air_C + params.dT_approach;
    % GEOPHIRES default injection temp is 70C, so we use it here.
    T_f_initial = 70;
    
    pathToScripts = 'GEOPHIRES\';
    
    % get reservoir impedance from reservoir
    result_reservoir = PorousReservoir(P_f_initial, T_f_initial, params.m_dot_IP, params);

    % Geophires only thinks it has one injection and production well
    % So for configurations with Field_IP_Multiplier > 1, we have to trick
    % it by using a reservoir impedance less than it actually is.
    % Also, reservoir impedance is calculated per IP pair, but geophires
    % puts the whole massflow through, so divide by the number of prod
    % wells
    reservoirImpedance = result_reservoir.ReservoirImpedance / Field_IP_Multiplier / numProdWells;
    
    reservoirFluidDensity = result_reservoir.AvgFluidDensity;
    m_dot_well = Field_IP_Multiplier * params.m_dot_IP;
    

    result_GEOPHIRES = RunGeophires(reservoirImpedance, reservoirFluidDensity, m_dot_well, numInjWells, numProdWells, params, pathToScripts);
    
    result.result_GEOPHIRES = result_GEOPHIRES;
    
    % set result
    % values returned by geophires are for entire system, not one IP
    result.Field_IP_Multiplier = Field_IP_Multiplier;
    result.reservoirImpedance = reservoirImpedance;
    result.W_pump_prod_IP = result_GEOPHIRES.W_pump_prod / Field_IP_Multiplier;
    result.T_prod_surface_C = result_GEOPHIRES.T_prod_surface_C;
    result.dP_pump = result_GEOPHIRES.dP_pump;
    result.W_net_IP = result_GEOPHIRES.W_net / Field_IP_Multiplier;
    
    % costs
    result.CapitalCost = result_GEOPHIRES.CapitalCost;
    if (result_GEOPHIRES.W_net < 0)
        result.SpecificCapitalCost_greenfield = NaN;
        result.LCOE_greenfield = NaN;
    else
        result.SpecificCapitalCost_greenfield = result_GEOPHIRES.CapitalCost / result_GEOPHIRES.W_net;
        result.LCOE_greenfield = result_GEOPHIRES.LCOE;
    end
    
    result.SpecificCapitalCost_brownfield = result.SpecificCapitalCost_greenfield;
    result.LCOE_brownfield = result.LCOE_greenfield;
end