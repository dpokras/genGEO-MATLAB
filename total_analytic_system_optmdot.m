function result = total_analytic_system_optmdot(params)

    tic;

    bestDependentVar = NaN;
    oldDependentVar = NaN;
    bestMdot = NaN;
    
    %Starting Massflow Guess
    params.m_dot_IP = 10;
        
    peaks = 0;
    d_m_dot = 10;
    d_m_dot_multiplier = 0.21; %0.25
    
    pointTowardsZero = false;
    changeMdotSign = false;
    
    % reversingFraction is the change past minimum to reverse
    reversingFraction = 0.03;
    
    while (peaks < 5)
        % Make sure mass flowrate is positive
        if (params.m_dot_IP <= 0)
            d_m_dot = -1 * d_m_dot * d_m_dot_multiplier;
            peaks = peaks + 1;
            params.m_dot_IP = abs(d_m_dot);
            disp(strcat(['mdot<=0; peak ' num2str(peaks)]));
        end
        
        % try at flowrate
        disp(strcat(['Trying a mass flowrate of ' num2str(params.m_dot_IP) ]));
        
        try
            % Run Simulation
            result = total_analytic_system(params);
      
            disp(strcat(['Production Temperature: ' num2str(result.T_prod_surface_C, 3) ...
                ' C; Power_IP: ' num2str(result.W_net_IP/1000, '%.1f') ...
                ' kWe; LCOE_brownfield: ' num2str(result.LCOE_brownfield*1e6, '%.1f') ...
                ' $/MWh; LCOE_greenfield: ' num2str(result.LCOE_greenfield*1e6, '%.1f') ...
                ' $/MWh.']));
            
            % Set Dependent Variable
            if (strcmp(params.optimizationMode, 'MinimizeLCOE_Brownfield'))
                dependentVar = result.LCOE_brownfield;
            elseif (strcmp(params.optimizationMode, 'MinimizeLCOE_Greenfield'))
                dependentVar = result.LCOE_greenfield;
            elseif (strcmp(params.optimizationMode, 'MaximizePower'))
                dependentVar = result.W_net_IP;
            else
                throw(MException('total_analytic_system_optmdot:NotImplemented','total_analytic_system_optmdot:OptimizationModeNotImplemented'));
            end
            
            % See if this mass flowrate is better than previous values
            if (isnan(dependentVar))
                %make sure it is trending towards mdot of zero
                pointTowardsZero = true;
            elseif (contains(params.optimizationMode, 'Minimize'))
                if (dependentVar < bestDependentVar || isnan(bestDependentVar))
                    bestDependentVar = dependentVar;
                    bestMdot = params.m_dot_IP;
                end
                if (peaks == 0)
                    % before first peak, exceed best val by
                    % reversingfraction
                    reversingThreshold = (bestDependentVar*(1+reversingFraction));
                    if (dependentVar > reversingThreshold)
                        changeMdotSign = true;
                        disp(strcat(['Minimize first crossing; peak ' num2str(peaks)]));
                    end
                else
                    % For subsequent peaks, be larger than previous value
                    if (dependentVar > oldDependentVar)
                        changeMdotSign = true;
                        disp(strcat(['Minimize crossing; peak ' num2str(peaks)]));
                    end
                end
            elseif (contains(params.optimizationMode, 'Maximize'))
                if (dependentVar > bestDependentVar || isnan(bestDependentVar))
                    bestDependentVar = dependentVar;
                    bestMdot = params.m_dot_IP;
                end
                if (peaks == 0)
                    % before first peak, exceed best val by
                    % reversingfraction
                    reversingThreshold = (bestDependentVar*(1-reversingFraction));
                    if (dependentVar < reversingThreshold)
                        changeMdotSign = true;
                        disp(strcat(['Maximize crossing; peak ' num2str(peaks)]));
                    end
                else
                    % For subsequent peaks, only be smaller than previous value
                    if (dependentVar < oldDependentVar)
                        changeMdotSign = true;
                        disp(strcat(['Maximize first crossing; peak ' num2str(peaks)]));
                    end
                end
            else
                throw(MException('total_analytic_system_optmdot:NotImplemented','total_analytic_system_optmdot:OptimizationModeNotImplemented'));
            end

            % Set residual
            if (isnan(dependentVar) || isnan(oldDependentVar) || d_m_dot == 0)
                residual = 0;
            else
                residual = abs(dependentVar - oldDependentVar)/abs(d_m_dot);
            end
            
            % Set current values to old values
            oldDependentVar = dependentVar;
               
        catch ME
            %catch some exceptions
            if (strcmp(ME.identifier,'MATLAB:Python:PyException') ...
                || strcmp(ME.identifier,'SWIG:RuntimeError') ...
                || strcmp(ME.identifier,'CapitalCost_SurfacePlant:NegativeTurbinePower') ...
                || strcmp(ME.identifier,'CapitalCost_SurfacePlant:NegativeBoilerHeat') ...
                || strcmp(ME.identifier,'CapitalCost_SurfacePlant:PositiveCondenserHeat') ...
                || strcmp(ME.identifier,'CapitalCost_SurfacePlant:PositiveCpgPumpPower') ...
                || strcmp(ME.identifier,'CapitalCost_SurfacePlant:PositiveProdPumpPower') ...
                || strcmp(ME.identifier,'total_analytic_system_CO2:TurbinePowerNegative') ...
                || strcmp(ME.identifier,'TotalAnalyticSystemWater:ExceedsMaxProductionPumpPressure') ...
                || strcmp(ME.identifier,'TotalAnalyticSystemWater:ExceedsMaxReservoirPressure') ...
                || strcmp(ME.identifier,'TotalAnalyticSystemCO2:ExceedsMaxReservoirPressure') ...
                || strcmp(ME.identifier,'TotalAnalyticSystemWater:ExceedsMaxTemperatureDecrease') ...
                || strcmp(ME.identifier,'ORC_Cycle_Tboil:T_in_NaN') ...
                || strcmp(ME.identifier,'ORC_Cycle_Tboil:T_boil_NaN') ...
                )
                %try again
                disp(strcat(['optimum_m_dot:Problem solving. ' ME.identifier ' Trying a different flowrate.']));
                
                % set these values to nan for the next loop.
                oldDependentVar = NaN;
                residual = 0;
                
                % make sure it goes towards zero
                pointTowardsZero = true;
           
            else
                % some other exception, stop simulation.
                disp(strcat(['optimum_m_dot:Exception: ' ME.message]));
                rethrow(ME)
            end
        end
        
        if (pointTowardsZero == true)
            pointTowardsZero = false;
            if (d_m_dot > 0)
                changeMdotSign = true;
                disp(strcat(['Pointing to zero.']));
            end
        end
        
        % Now adjust massflow
        if (changeMdotSign == true)
            changeMdotSign = false;
            d_m_dot = -1 * d_m_dot * d_m_dot_multiplier;
            peaks = peaks + 1;
            disp(strcat(['Crossed peak ' num2str(peaks) ', Changing direction.']));
        end
        
        % set new massflowrate
        params.m_dot_IP = params.m_dot_IP + d_m_dot; 
            
    end
    
    
    
    
    
    
    % The code below reruns the simulation at the best mdot and returns the
    % results.
    
    
    
    % if no solution was ever found, return all defaults.
    % otherwise, return results for either the best LCOE or mdot.
    if (isnan(bestMdot)==true)
        % Set all output variables to their default (0 or NaN)
        result.dP_turbine = NaN;
        result.T_prod_surface_C = NaN;
        result.W_net_IP = 0;
        result.Field_IP_Multiplier = 0;
        result.Q_fluid_IP = 0;
        result.N_IP_Multiplier = 0;
        result.State_injection = [];
        result.State_reservoir = [];
        result.State_production = [];
        result.CapitalCost.CostSurfacePlant.C_T_G = NaN;
        result.CapitalCost.CostSurfacePlant.C_pump_orc = NaN;
        result.CapitalCost.CostSurfacePlant.C_coolingTowers = NaN;
        result.CapitalCost.CostSurfacePlant.C_heatExchanger = NaN;
        result.CapitalCost.CostSurfacePlant.C_recuperator = NaN;
        result.CapitalCost.CostSurfacePlant.C_pump_prod = NaN;
        result.CapitalCost.CostSurfacePlant.C_pump_inj = NaN;
        result.CapitalCost.CostSurfacePlant.C_plant_otherEquipment = NaN;
        result.CapitalCost.CostSurfacePlant.C_plant_installation = NaN;
        result.CapitalCost.CostSurfacePlant.C_plant_indirectContingency = NaN;
        result.CapitalCost.C_surfacePlant = NaN;
        result.CapitalCost.C_gatheringSystem = NaN;
        result.CapitalCost.C_wells_production = NaN;
        result.CapitalCost.C_wells_injection = NaN;
        result.CapitalCost.C_wellfield = NaN;
        result.CapitalCost.C_exploration = NaN;
        result.CapitalCost.C_stimulation = NaN;
        result.CapitalCost.C_brownfield = NaN;
        result.CapitalCost.C_greenfield = NaN;
        result.SpecificCapitalCost_brownfield = NaN;
        result.SpecificCapitalCost_greenfield = NaN;
        result.LCOE_brownfield = NaN;
        result.LCOE_greenfield = NaN;
        
        result.m_dot_IP = 0;
        result.residual = NaN;
  
        disp('Found a min LCOE and Power of zero at zero flow.');
    else
        %reevaluate system at best flowrate
        if (strcmp(params.optimizationMode, 'MinimizeLCOE_Brownfield'))
            disp(strcat(['Found a min LCOE_Brownfield of ' num2str(bestDependentVar*1e6, '%.0f') ' $/MWh at a mass flowrate of ' num2str(bestMdot, '%.1f') ' kg/s.']));
        elseif (strcmp(params.optimizationMode, 'MinimizeLCOE_Greenfield'))
            disp(strcat(['Found a min LCOE_Greenfield of ' num2str(bestDependentVar*1e6, '%.0f') ' $/MWh at a mass flowrate of ' num2str(bestMdot, '%.1f') ' kg/s.']));
        elseif (strcmp(params.optimizationMode, 'MaximizePower'))            
            disp(strcat(['Found a max power of ' num2str(bestDependentVar/1e3, '%.1f') ' kWe at a mass flowrate of ' num2str(bestMdot, '%.1f') ' kg/s.']));
        else
            disp('mode not found');
        end
        
        params.m_dot_IP = bestMdot;
        
        disp(strcat(['Final simulation at a mass flowrate of ' num2str(params.m_dot_IP, '%.4f') ' kg/s.']));
        result = total_analytic_system(params);
        result.m_dot_IP = params.m_dot_IP;
        result.residual = residual;
    end
    
    result.secondsToSolve = toc;
    
end