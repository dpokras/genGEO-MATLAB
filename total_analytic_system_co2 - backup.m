function result = total_analytic_system_co2(params)
    
    if (strcmp(params.fluid,'CO2') == false)
        throw(MException('total_analytic_system_CO2:Wrong_Fluid','Wrong Fluid--Not CO2'));
    end
    
    % Get condensation pressure
    T_condensation = params.T_surface_air_C + params.dT_approach;
    P_condensation = CoolProp.PropsSI('P', 'T', T_condensation+273.15, 'Q', 0, params.fluid) + 50e3;
    % Start without pumping
    dP_pump = 0;
    
    P_reservoir = params.depth * 1000 * 9.81;
    P_reservoir_max = params.depth * 2500 * 9.81;
    
    % If a porous system, loop until the downhole pressure is near
    % hydrostatic
    dP_downhole_threshold = 1e3; %1 kPa
    dP_downhole = NaN;
    dP_loops = 1;
    dP_Solver = Solver;
    %dP_Solver.showPlot = true;
    % Dependent Variable for this solver is dP_pump
    % Independent Variable for this solver is dP_downhole
    while (isnan(dP_downhole) || abs(dP_downhole) >= dP_downhole_threshold)
    
        %Find Injection Conditions
        % if dP_pump is less than zero, throttle downhole to avoid twophase
        % at the surface
        P_pump_inlet = P_condensation;
        if (dP_pump >= 0)
            P_pump_outlet = P_pump_inlet + dP_pump;
        else
            P_pump_outlet = P_pump_inlet;
        end
        T_pump_inlet = T_condensation;
        h_pump_inlet = CoolProp.PropsSI('HMASS', 'P', P_pump_inlet, 'T', T_pump_inlet+273.15, params.fluid);
        s_pump_inlet = CoolProp.PropsSI('SMASS', 'P', P_pump_inlet, 'T', T_pump_inlet+273.15, params.fluid);
    
        % if dP_pump positive, it's a pump, if negative it's a throttle.
        if (dP_pump > 0)
            h_pump_outletS = CoolProp.PropsSI('HMASS', 'P', P_pump_outlet, 'S', s_pump_inlet, params.fluid);
            h_pump_outlet = h_pump_inlet + (h_pump_outletS - h_pump_inlet) / params.eta_cpg_pump;
        else
            h_pump_outlet = h_pump_inlet;
        end
    
        W_pump = -1 * params.m_dot * (h_pump_outlet - h_pump_inlet);
    
        P_inj_surface = P_pump_outlet;
        T_inj_surface = CoolProp.PropsSI('T', 'P', P_pump_outlet, 'HMASS', h_pump_outlet, params.fluid) - 273.15;
    
        %injection well
    
        %%%%%%%
        %%%%%%%
        %%%%%%%
    
    
    
        %%%%%%% Add a for loop here to iterate through different distances
        %%%%%%% between side streams. Collect the properties for different
        %%%%%%% depths and calculate the whole loop n_stream times.
    
        depths = params.depth - params.stream_spacing*((params.n_streams-1):-1:0);

        temp_at_depth = params.T_surface_rock_C + params.dT_dz * depths; %C

        %%%%%%%
        %%%%%%%
        %%%%%%%
        disp('simulating injection well')
        result_injWell = semi_analytic_well(P_inj_surface, T_inj_surface, params.T_surface_rock_C, -1*depths, 0, params.m_dot, params.well_radius, params);
        
        %after injection well, throttle if dP_pump negative
        %dP_pump isn't often negative, usually in shallow, low flowrate
        %cases
        if (dP_pump < 0)
            P_inj_downhole = result_injWell.EndPressure + dP_pump;
        else
            P_inj_downhole = result_injWell.EndPressure;
        end

        h_inj_downhole = result_injWell.EndEnthalpy;
        T_inj_downhole = CoolProp.PropsSI('T', 'P', P_inj_downhole, 'HMASS', h_inj_downhole, params.fluid) - 273.15;

        %reservoir
        disp('simulating reservior')

        if (strcmp(params.system,'Porous'))
            result_reservoir = PorousReservoir(P_inj_downhole, T_inj_downhole, params.m_dot, params);
            % find downhole pressure difference (negative means
            % overpressure
            dP_downhole = P_reservoir - result_reservoir.EndPressure;
        else
            result_reservoir = semi_analytic_well(P_inj_downhole, T_inj_downhole, temp_at_depth, 0, params.res_length, params.m_dot/params.n_streams, params.side_stream_radius, params);
            % No need to loop if using closed reservoir
            dP_downhole = 0;
        end

        

    
        % calculate what next
        % independentVariable_zeroGuess = Solver(independentVariable, dependentVariable);
        dP_pump = dP_Solver.AddDataAndEstimate(dP_pump, dP_downhole);
        % dP_pump can be nan if only one trial has been done
        if (isnan(dP_pump))
            dP_pump = 0.5 * dP_downhole;
        end

        % Warn against excessive loops
        if (dP_loops > 10)
            disp(strcat(['Warning: total_analytic_system_water:dP1_loops is large: ' num2str(dP_loops)]));
        end
        dP_loops = dP_loops + 1;
    end
    
    % if reservoir pressure greater than allowable, throw error.
    if (result_reservoir.EndPressure(end) >= P_reservoir_max)
        throw(MException('TotalAnalyticSystemCO2:ExceedsMaxReservoirPressure',strcat(['Exceeds Max Reservoir Pressure of ' num2str(P_reservoir_max/1e6,3) ' MPa!'])));
    end

    %production well and surface result

%     running_depth = params.depth - 100*(iter-1);
%     temp_at_depth = params.T_surface_rock_C + params.dT_dz * running_depth; %C
    disp('simulating production well')
    result_prodWell = semi_analytic_well(result_reservoir.EndPressure, result_reservoir.EndTemp, temp_at_depth, depths, 0, params.m_dot, params.well_radius, params);
    
    %%%% Take the average temperature arriving at the production well
        %%% *** only if the stream flowrates are equal***
        %%% *** otherwise average proportionally ***
    T_prod_surface = mean([result_prodWell.EndTemp]);
    %%%% Take the lowest pressure of the separate side streams
    P_prod_surface = result_prodWell.EndPressure(1);    
    
    %%%%%%
    %%%%%%
    %%%%%%
    %%%%%% Adjust calculation of friction in pipes by seperating the
    %%%%%% sections for side streams and main well
    %%%%%%
    %%%%%%
    %%%%%%
    
    %%%%%%
    %%%%%% VVVVVVV Where is this coming from??? VVVVVVV
    %%%%%%
    
    % pressure losses between production well and surface plant
    P_surfacePlant = Surface_plant.plant_dP(P_prod_surface, T_prod_surface, params);
    T_prod_surface
    %Calculate turbine work
    result_turbine = Surface_plant.TurbWorkfunc(P_surfacePlant, T_prod_surface, params);
    W_turbine = result_turbine.W_turbine;

    %Calculate condenser work
    [W_cond_HX, Q_cond_HX] = Surface_plant.CondWorkfunc(result_turbine.P, result_turbine.T, params);
%     W_cond_HX = 0;
%     Q_cond_HX = 0;
    % Calculate net work
    W_net = Surface_plant.NetWorkfunc(W_turbine, W_cond_HX, W_pump);

    %%%%%%
    %%%%%% Changing P_condenstion to be an iteration of P until liquid
    %%%%%% condensation.
    %%%%%%
    
    %end

    Q_fluid = (result_injWell.Heat + result_reservoir.Heat + result_prodWell.Heat);
    
    %Calculate Capital Cost
    result_capitalCost = CapitalCost_CPG(W_turbine, 0, Q_cond_HX, W_pump, 0, params);
    result_capitalCost.CostSurfacePlant.C_pump_orc = 0;
    result_capitalCost.CostSurfacePlant.C_heatExchanger = 0;
    result_capitalCost.CostSurfacePlant.C_recuperator = 0;
    result_capitalCost.CostSurfacePlant.C_pump_prod = 0;
    
    %Calculate LCOE ($/Weh)
    %Calculate Specific Capital Cost ($/We)
    result_brownfield = LCOE_Simple(result_capitalCost.C_brownfield, W_net, params);
    result_greenfield = LCOE_Simple(result_capitalCost.C_greenfield, W_net, params);
    
    
    % set result
    % set result_wells
    result.injection = result_injWell;
    result.reservoir = result_reservoir;
    result.productionLower = [];
    result.productionUpper = result_prodWell;
    
    % set result_plant
    result.Q_preheater = 0;
    result.Q_boiler = 0;
    result.W_turbine = W_turbine;
    result.W_cond_HX = W_cond_HX;
    result.Q_recuperator = 0;
    result.Q_desuperheater = 0;
    result.Q_cond_HX = Q_cond_HX;
    result.W_pump_orc = 0;
    result.W_cooler_orc = 0;
    result.W_condenser_orc = 0;
    result.W_pump_prod = 0;
    result.W_pump = W_pump;
    result.W_net = W_net;
    
    % set result_fluid
    result.maxPump_dP = 0;
    result.maxSurface_dT = 0;
    result.temp_at_depth = temp_at_depth;
    result.P_reservoir = P_reservoir;
    result.P_reservoir_max = P_reservoir_max;
    result.dT_inj = 0;
    result.dP_downhole = dP_downhole;
    result.dP_pump = dP_pump;
%     result.dP_turbine = dP_turbine;
%     result.dP_surfacePipes = dP_surfacePipes;
    result.dP_downhole_loops = dP_loops;
    
    result.T_prod_surface = T_prod_surface;
    result.W_net = W_net;
    result.Q_fluid = Q_fluid;
    result.max_speed = params.m_dot/(pi*(params.well_radius)^2)/result_prodWell.EndDensity(1);
%     result.Field_Multiplier = Field_Multiplier;
    
    % costs
    result.CapitalCost = result_capitalCost;
    result.SpecificCapitalCost_brownfield = result_brownfield.SpecificCapitalCost;
    result.SpecificCapitalCost_greenfield = result_greenfield.SpecificCapitalCost;
    result.LCOE_brownfield = result_brownfield.LCOE;
    result.LCOE_greenfield = result_greenfield.LCOE;
end