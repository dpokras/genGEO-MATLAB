function result = total_analytic_system_water(params)

    %Check that simulation is for water
    if (strcmp(params.fluid,'Water') == false)
        throw(MException('total_analytic_system_water:Wrong_Fluid','Wrong Fluid--Not Water'));
    end
    
    %Prod/Inj Well Multipliers
    if (strcmp(params.wellFieldType, 'Doublet'))
        injectionWellMultiplier = 1;
        productionWellMultiplier = 1;
        Field_IP_Multiplier = 1;
    elseif (strcmp(params.wellFieldType, '5spot_SharedNeighbor'))
        injectionWellMultiplier = 4;
        productionWellMultiplier = 4;
        Field_IP_Multiplier = 4;
    elseif (strcmp(params.wellFieldType, '5spot'))
        injectionWellMultiplier = 4;
        productionWellMultiplier = 1;
        Field_IP_Multiplier = 4;
    elseif (strcmp(params.wellFieldType, '5spot_ManyN'))
        injectionWellMultiplier = 4;
        if (N_5spot == 1)
            productionWellMultiplier = 1;
        else
            productionWellMultiplier = 4;
        end
        Field_IP_Multiplier = 4 * N_5spot^2;
    elseif (strcmp(params.wellFieldType, 'Tungsten'))
        injectionWellMultiplier = 1;
        productionWellMultiplier = 1;
        Field_IP_Multiplier = 4;
    else
        throw(MException('total_analytic_system_water:BadWellfieldType','Bad Wellfield Type'));
    end
    
   
    temp_at_depth = params.T_surface_rock_C + params.dT_dz * params.depth; %C
    P_reservoir = params.depth * 1000 * 9.81;
    P_reservoir_max = params.depth * 2500 * 9.81;
    
    P_system_min = CoolProp('P', 'T', temp_at_depth+273.15, 'Q', 0, 'Water') + 1e5;

    %Set Injection Conditions
    % Guess two numbers to start
    P_inj_surface_initial = P_system_min + 1e6; %Guess 1MPa above flash point first
    T_inj_surface_initial_C = 60;
    h_inj_surface_initial = CoolProp('HMASS', 'P', P_inj_surface_initial, 'T', T_inj_surface_initial_C+273.15, params.fluid);

    % calculate surface pipe loss friction factor
    ff = FrictionFactor(params.well_radius, P_inj_surface_initial, h_inj_surface_initial, params.m_dot_IP, params);
    
    % loop until water injection temperature is determined
    % This loop solves for T_inj_surface_C (independent variable) to make
    % dT_inj (dependent variable) zero.
    T_inj_surface_C = T_inj_surface_initial_C;
    dT_inj = NaN;      
    dT_loops = 1;
    % function used to drive dT_inj to zero.
    dT_Solver = Solver;
    %dT_Solver.showPlot = true;
    while (isnan(dT_inj) || abs(dT_inj) >=0.5)
        
        % Find necessary injection pressure
        % This loop solves for P_inj_surface (independent variable) to make
        % dP_downhole (dependent variable) zero.
        P_inj_surface = P_inj_surface_initial;
        dP_downhole = NaN;
        dP1_loops = 1;
        dP1_Solver = Solver;
        %dP1_Solver.showPlot = true;
        while ( isnan(dP_downhole) || abs(dP_downhole) > 10e3 )
            %injection well
            result_injWell = semi_analytic_well(P_inj_surface, T_inj_surface_C, params.T_surface_rock_C, -1*params.depth, 0, params.m_dot_IP*injectionWellMultiplier, params.well_radius, params);
            
            %reservoir well
            if (strcmp(params.system,'Porous'))
                result_reservoir = PorousReservoir(result_injWell.EndPressure, result_injWell.EndTemp, params.m_dot_IP, params);
            elseif (strcmp(params.system,'Conduction1'))
                result_reservoir = semi_analytic_well(result_injWell.EndPressure, result_injWell.EndTemp, temp_at_depth, 0, params.res_length, params.m_dot_IP, params.well_radius, params);
            elseif (strcmp(params.system,'Conduction2'))
                result_reservoir = semi_analytic_well(result_injWell.EndPressure, result_injWell.EndTemp, temp_at_depth, 0, params.res_length, params.m_dot_IP/2, params.well_radius, params);
            elseif (strcmp(params.system,'Conduction4'))
                result_reservoir = semi_analytic_well(result_injWell.EndPressure, result_injWell.EndTemp, temp_at_depth, 0, params.res_length, params.m_dot_IP/4, params.well_radius, params);
            elseif (strcmp(params.system,'Conduction8'))
                result_reservoir = semi_analytic_well(result_injWell.EndPressure, result_injWell.EndTemp, temp_at_depth, 0, params.res_length, params.m_dot_IP/8, params.well_radius, params);
            else
                throw(MException('totalSystem:NotImplemented','total_analytic_system_water:Not Implemented'));
            end

            %if already at P_system_min, stop looping
            if (P_inj_surface == P_system_min)
                break;
            end
            
            % find downhole pressure difference (negative means overpressure)
            dP_downhole = P_reservoir - result_reservoir.EndPressure;
            P_inj_surface = dP1_Solver.AddDataAndEstimate(P_inj_surface, dP_downhole);
            if (isnan(P_inj_surface))
                P_inj_surface = P_inj_surface_initial + dP_downhole;
            end
            
            % Warn against excessive loops
            if (dP1_loops > 10)
                disp(strcat(['Warning: total_analytic_system_water:dP1_loops is large: ' num2str(dP1_loops)]));
            end
            dP1_loops = dP1_loops + 1;
            
            % Set Limit
            if (P_inj_surface < P_system_min)
                % can't be below this temp or fluid will flash
                P_inj_surface = P_system_min;
                %disp('total_analytic_system_water:P_inj_surface was below flash pressure! Correcting.');
            end
        end
        
        % if reservoir pressure greater than allowable, throw error.
        if (result_reservoir.EndPressure >= P_reservoir_max)
            throw(MException('TotalAnalyticSystemWater:ExceedsMaxReservoirPressure',strcat(['Exceeds Max Reservoir Pressure of ' num2str(P_reservoir_max/1e6,3) ' MPa!'])));
        end
        
        %production well (first of two wells)
        result_prodWellLower = semi_analytic_well(result_reservoir.EndPressure, result_reservoir.EndTemp, temp_at_depth, (params.depth-params.PumpDepth), 0, ...
            params.m_dot_IP*productionWellMultiplier, params.well_radius, params);
        
              
        % Now pumping and second well
        % This loop solves for dP_pump (independent variable) to make
        % dP_surface (dependent variable) zero.
        dP_pump = 0;
        dP_surface = NaN;
        dP2_loops = 1;
        dP2_Solver = Solver;
        %dP2_Solver.showPlot = true;
        while (isnan(dP_surface) || abs(dP_surface) > 100)
            try
                if (dP_pump > params.maxPump_dP)
                    dP_pump = params.maxPump_dP;
                end
                P_prod_pump_out = result_prodWellLower.EndPressure + dP_pump;
                T_prod_pump_out = result_prodWellLower.EndTemp;
                temp_at_pump_depth = params.T_surface_rock_C + params.dT_dz * params.PumpDepth;

                result_prodWellUpper = semi_analytic_well(P_prod_pump_out, T_prod_pump_out, temp_at_pump_depth, params.PumpDepth, 0, ...
                    params.m_dot_IP*productionWellMultiplier, params.well_radius, params);
                
                % calculate surface losses
                if (params.hasSurfaceGatheringSystem == true)
                    dP_surfacePipes = ff * params.res_length / (2*params.well_radius)^5 * 8 * params.m_dot_IP^2 / result_prodWellUpper.EndDensity / pi^2;
                else
                    dP_surfacePipes = 0;
                end

                dP_surface = (result_prodWellUpper.EndPressure - dP_surfacePipes - P_inj_surface);
                
                if (dP_pump == params.maxPump_dP)
                    break;
                end
                
                % No pumping is needed in this system
                if (dP_pump == 0 && dP_surface >= 0)
                    break;
                end
                
                dP_pump = dP2_Solver.AddDataAndEstimate(dP_pump, dP_surface);
                if (isnan(dP_pump))
                    dP_pump = -1 * dP_surface;
                end
                
                % Pump can't be less than zero
                if (dP_pump < 0)
                    %disp('total_analytic_system_water:dP_pump was less than zero! Correcting.');
                    dP_pump = 0;
                end
                
                % Warn against excessive loops
                if (dP2_loops > 10)
                    disp(strcat(['Warning: total_analytic_system_water:dP2_loops is large: ' num2str(dP2_loops)]));
                end
                dP2_loops = dP2_loops + 1;

            catch ME
                % Only catch problems of flashing fluid
               if (strcmp(ME.identifier,'SemiAnalytic:BelowSaturationPressure'))
                   % slowly increase pressure
                   dP_pump = dP_pump + 1e5;
               else
                    rethrow(ME)
               end
            end 
        end
        
        % if pump pressure greater than allowable, throw error
        if (dP_pump >= params.maxPump_dP)
            throw(MException('TotalAnalyticSystemWater:ExceedsMaxProductionPumpPressure',strcat(['Exceeds Max Pump Pressure of ' num2str(params.maxPump_dP/1e6,3) ' MPa!'])));
        end
        
        % Calculate pump power
        if (dP_pump > 0)
            h_prod_pump_out = CoolProp('HMASS', 'P', P_prod_pump_out, 'T', T_prod_pump_out+273.15, params.fluid);
            W_pump_prod_IP = params.m_dot_IP * (result_prodWellLower.EndEnthalpy - h_prod_pump_out) / params.eta_pump;
        else
            % just a throttle, no power
            W_pump_prod_IP = 0;
        end
        
        % Calculate injection temp
        T_prod_surface_C = result_prodWellUpper.EndTemp;
        result_ORC = ORC_Cycle(T_prod_surface_C, params);
        
        dT_inj = T_inj_surface_C - result_ORC.T_out_C;
        
        %Solver
        T_inj_surface_C = dT_Solver.AddDataAndEstimate(T_inj_surface_C, dT_inj);
        if (isnan(T_inj_surface_C))
            T_inj_surface_C = result_ORC.T_out_C;
        end
        % add lower bounds
        if (T_inj_surface_C < 1)
            T_inj_surface_C = 1;
        end
        % add upper bounds that make sense
        if (T_inj_surface_C > T_prod_surface_C && T_inj_surface_C > 50)
            T_inj_surface_C = T_prod_surface_C;
        end
        
        % if excessive loops, alert (we don't want an endless looping
        % situation)
        if (dT_loops > 10)
            disp(strcat(['Warning: total_analytic_system_water:dT_loops is large: ' num2str(dT_loops)]));
        end
        dT_loops = dT_loops + 1;
        %disp(strcat(['SolverStatus_IVs: P_inj_surface=' num2str(P_inj_surface) ', dP_pump=' num2str(dP_pump) ', T_inj_surface_C=' num2str(T_inj_surface_C)]));
        %disp(strcat(['SolverStatus_DVs: dP_downhole=' num2str(dP_downhole) ', dP_surface=' num2str(dP_surface) ', dT_inj=' num2str(dT_inj)]));
    end
    
    % Check injection temperature
    if (strcmp(params.silicaPrecip,'PreventSilica'))
        maxSurface_dT = 89; % To prevent silica precip...from dipippo 1985
        if (T_prod_surface_C - T_inj_surface_C > params.maxSurface_dT)
            throw(MException('TotalAnalyticSystemWater:ExceedsMaxTemperatureDecrease',strcat(['Exceeds Max Temp Decrease of ' num2str(maxSurface_dT,3) ' C to prevent silica precipitation!'])));
        end
    elseif (strcmp(params.silicaPrecip,'IgnoreSilica'))
        maxSurface_dT = Inf;
        % OK
    else
        throw(MException('total_analytic_system_water:UnknownSilica','Silica Option Unknown'));
    end
    
    
    % Power Plant Models
    if (strcmp(params.orcModel,'Simulation'))
        result_ORC = ORC_Cycle(T_prod_surface_C, params);
    elseif (strcmp(params.orcModel,'Geophires'))
        result_ORC = ORC_Cycle_GEOPHIRES(T_prod_surface_C, params);
        result_ORC.q_preheater = 0;
        result_ORC.q_recuperator = 0;
        result_ORC.q_desuperheater = 0;
        result_ORC.w_pump = 0;
        result_ORC.w_cooler = 0;
        result_ORC.w_condenser = 0;
        result_ORC.dT_range_CT = 7;
        result_ORC.dT_LMTD_preheater = 7;
        result_ORC.dT_LMTD_boiler = 7;
        result_ORC.dT_LMTD_recuperator = NaN;
    else
        throw(MException('total_analytic_system_water:UnknownOrcModel','Unknown ORC Model'));
    end
    
    Q_preheater_IP = params.m_dot_IP * result_ORC.q_preheater;
    Q_boiler_IP = params.m_dot_IP * result_ORC.q_boiler;
    W_turbine_IP = params.m_dot_IP * result_ORC.w_turbine;
    Q_recuperator_IP = params.m_dot_IP * result_ORC.q_recuperator;
    Q_desuperheater_IP = params.m_dot_IP * result_ORC.q_desuperheater;
    Q_condenser_IP = params.m_dot_IP * result_ORC.q_condenser;
    W_pump_orc_IP = params.m_dot_IP * result_ORC.w_pump;
    W_cooler_orc_IP = params.m_dot_IP * result_ORC.w_cooler;
    W_condenser_orc_IP = params.m_dot_IP * result_ORC.w_condenser;
    
    W_net_IP = W_turbine_IP + W_pump_orc_IP + W_cooler_orc_IP + W_condenser_orc_IP + W_pump_prod_IP;
    Q_fluid_IP = (Q_preheater_IP + Q_boiler_IP);
    
    % Caculate total for system
    Q_preheater_total = Field_IP_Multiplier * Q_preheater_IP;
    Q_boiler_total = Field_IP_Multiplier * Q_boiler_IP;
    W_turbine_total = Field_IP_Multiplier * W_turbine_IP;
    Q_recuperator_total = Field_IP_Multiplier * Q_recuperator_IP;
    Q_desuperheater_total = Field_IP_Multiplier * Q_desuperheater_IP;
    Q_condenser_total = Field_IP_Multiplier * Q_condenser_IP;
    W_pump_orc_total = Field_IP_Multiplier * W_pump_orc_IP;
    W_pump_prod_total = Field_IP_Multiplier * W_pump_prod_IP;
    W_net_total = Field_IP_Multiplier * W_net_IP;
    
    % Find Capital Cost
    result_capitalCost = CapitalCost_Water(Q_preheater_total, Q_boiler_total, W_turbine_total, Q_recuperator_total, Q_desuperheater_total, ...
        Q_condenser_total, W_pump_orc_total, W_pump_prod_total, ...
        result_ORC.dT_range_CT, result_ORC.dT_LMTD_preheater, result_ORC.dT_LMTD_boiler, result_ORC.dT_LMTD_recuperator, params);
    result_capitalCost.CostSurfacePlant.C_pump_inj = 0;
    
    %Calculate LCOE ($/We-h)
    %Calculate Specific Capital Cost ($/We)
    result_brownfield = LCOE_Simple(result_capitalCost.C_brownfield, W_net_total, params);
    result_greenfield = LCOE_Simple(result_capitalCost.C_greenfield, W_net_total, params);
    
    % set result_wells
    result.injection = result_injWell;
    result.reservoir = result_reservoir;
    result.productionLower = result_prodWellLower;
    result.productionUpper = result_prodWellUpper;
    
    % set result_plant
    result.result_ORC = result_ORC;
    result.Q_preheater_IP = Q_preheater_IP;
    result.Q_boiler_IP = Q_boiler_IP;
    result.W_turbine_IP = W_turbine_IP;
    result.Q_recuperator_IP = Q_recuperator_IP;
    result.Q_desuperheater_IP = Q_desuperheater_IP;
    result.Q_condenser_IP = Q_condenser_IP;
    result.W_pump_orc_IP = W_pump_orc_IP;
    result.W_cooler_orc_IP = W_cooler_orc_IP;
    result.W_condenser_orc_IP = W_condenser_orc_IP;
    result.W_pump_prod_IP = W_pump_prod_IP;
    result.W_pump_inj_IP = 0;
    result.W_net_IP = W_net_IP;
    
    % set result_fluid
    result.maxSurface_dT = maxSurface_dT;
    result.temp_at_depth = temp_at_depth;
    result.P_reservoir = P_reservoir;
    result.P_reservoir_max = P_reservoir_max;
    result.P_system_min = P_system_min;
    result.dT_inj = dT_inj;
    result.dP_downhole = dP_downhole;
    result.dP_pump = dP_pump;
    result.dP_turbine = 0;
    result.dP_surfacePipes = dP_surfacePipes;
    result.dP1_downhole_loops = dP1_loops;
    result.dP2_downhole_loops = dP2_loops;
    result.dT_loops = dT_loops;
    
    result.T_prod_surface_C = T_prod_surface_C;
    result.T_inj_surface_C = T_inj_surface_C;
    result.W_net_IP = W_net_IP;
    result.Q_fluid_IP = Q_fluid_IP;
    result.Field_IP_Multiplier = Field_IP_Multiplier;
    
    % costs
    result.CapitalCost = result_capitalCost;
    result.LCOE_brownfield = result_brownfield.LCOE;
    result.SpecificCapitalCost_brownfield = result_brownfield.SpecificCapitalCost;
    result.LCOE_greenfield = result_greenfield.LCOE;
    result.SpecificCapitalCost_greenfield = result_greenfield.SpecificCapitalCost;
end