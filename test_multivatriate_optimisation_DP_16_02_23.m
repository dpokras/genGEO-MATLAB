clear;
close all;clc
% Purpose of sheet is to generate a plot for the different system
% configurations (1 Lateral, 2 Laterals,...) for a variety of reservoir
% depths
%% Define System Parameters
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
params.wellFieldType = 'Doublet'; 
params.time_years = 30;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';

% Select Fluid Type
fluid_All = ["Water","CO2"];
params.fluid = fluid_All(2);
%% Base Case
params.dT_dz = 0.035;
params.wellCostType = 'Baseline';

prob = optimproblem;

% Optimisation variables
params.m_dot_IP = optimvar('mdots', 'LowerBound', 0);
params.res_length = optimvar('res_length', 'LowerBound',500, 'UpperBound',5000);
params.depth = optimvar('depth', 'LowerBound', 1000, 'UpperBound', 8000);
% params.well_radius = optimvar('well_radius', 'LowerBound', 0.1, 'UpperBound', 0.5);
params.n_streams = optimvar('n_streams', 'LowerBound', 1, 'UpperBound',20);

% Run Simulation
result = total_analytic_system(params);


% Contraints

constr = result.W_net_IP == 1300;
prob.Constraints.powerconstr = constr;

% Objective function

cost_g = result.LCOE_greenfield;

prob.Objective = cost_g;

showproblem(prob)

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
%             
%             % See if this mass flowrate is better than previous values
%             if (isnan(dependentVar))
%                 %make sure it is trending towards mdot of zero
%                 pointTowardsZero = true;
%             elseif (contains(params.optimizationMode, 'Minimize'))
%                 if (dependentVar < bestDependentVar || isnan(bestDependentVar))
%                     bestDependentVar = dependentVar;
%                     bestMdot = params.m_dot_IP;
%                 end
%                 if (peaks == 0)
%                     % before first peak, exceed best val by
%                     % reversingfraction
%                     reversingThreshold = (bestDependentVar*(1+reversingFraction));
%                     if (dependentVar > reversingThreshold)
%                         changeMdotSign = true;
%                         disp(strcat(['Minimize first crossing; peak ' num2str(peaks)]));
%                     end
%                 else
%                     % For subsequent peaks, be larger than previous value
%                     if (dependentVar > oldDependentVar)
%                         changeMdotSign = true;
%                         disp(strcat(['Minimize crossing; peak ' num2str(peaks)]));
%                     end
%                 end
%             elseif (contains(params.optimizationMode, 'Maximize'))
%                 if (dependentVar > bestDependentVar || isnan(bestDependentVar))
%                     bestDependentVar = dependentVar;
%                     bestMdot = params.m_dot_IP;
%                 end
%                 if (peaks == 0)
%                     % before first peak, exceed best val by
%                     % reversingfraction
%                     reversingThreshold = (bestDependentVar*(1-reversingFraction));
%                     if (dependentVar < reversingThreshold)
%                         changeMdotSign = true;
%                         disp(strcat(['Maximize crossing; peak ' num2str(peaks)]));
%                     end
%                 else
%                     % For subsequent peaks, only be smaller than previous value
%                     if (dependentVar < oldDependentVar)
%                         changeMdotSign = true;
%                         disp(strcat(['Maximize first crossing; peak ' num2str(peaks)]));
%                     end
%                 end
%             else
%                 throw(MException('total_analytic_system_optmdot:NotImplemented','total_analytic_system_optmdot:OptimizationModeNotImplemented'));
%             end
% 
%             % Set residual
%             if (isnan(dependentVar) || isnan(oldDependentVar) || d_m_dot == 0)
%                 residual = 0;
%             else
%                 residual = abs(dependentVar - oldDependentVar)/abs(d_m_dot);
%             end
%             
%             % Set current values to old values
%             oldDependentVar = dependentVar;
%                
%         catch ME
%             %catch some exceptions
%             if (strcmp(ME.identifier,'MATLAB:Python:PyException') ...
%                 || strcmp(ME.identifier,'SWIG:RuntimeError') ...
%                 || strcmp(ME.identifier,'CapitalCost_SurfacePlant:NegativeTurbinePower') ...
%                 || strcmp(ME.identifier,'CapitalCost_SurfacePlant:NegativeBoilerHeat') ...
%                 || strcmp(ME.identifier,'CapitalCost_SurfacePlant:PositiveCondenserHeat') ...
%                 || strcmp(ME.identifier,'CapitalCost_SurfacePlant:PositiveCpgPumpPower') ...
%                 || strcmp(ME.identifier,'CapitalCost_SurfacePlant:PositiveProdPumpPower') ...
%                 || strcmp(ME.identifier,'total_analytic_system_CO2:TurbinePowerNegative') ...
%                 || strcmp(ME.identifier,'TotalAnalyticSystemWater:ExceedsMaxProductionPumpPressure') ...
%                 || strcmp(ME.identifier,'TotalAnalyticSystemWater:ExceedsMaxReservoirPressure') ...
%                 || strcmp(ME.identifier,'TotalAnalyticSystemCO2:ExceedsMaxReservoirPressure') ...
%                 || strcmp(ME.identifier,'TotalAnalyticSystemWater:ExceedsMaxTemperatureDecrease') ...
%                 || strcmp(ME.identifier,'ORC_Cycle_Tboil:T_in_NaN') ...
%                 || strcmp(ME.identifier,'ORC_Cycle_Tboil:T_boil_NaN') ...
%                 )
%                 %try again
%                 disp(strcat(['optimum_m_dot:Problem solving. ' ME.identifier ' Trying a different flowrate.']));
%                 
%                 % set these values to nan for the next loop.
%                 oldDependentVar = NaN;
%                 residual = 0;
%                 
%                 % make sure it goes towards zero
%                 pointTowardsZero = true;
%            
%             else
%                 % some other exception, stop simulation.
%                 disp(strcat(['optimum_m_dot:Exception: ' ME.message]));
%                 rethrow(ME)
%             end
%         end
%         
%         if (pointTowardsZero == true)
%             pointTowardsZero = false;
%             if (d_m_dot > 0)
%                 changeMdotSign = true;
%                 disp(strcat(['Pointing to zero.']));
%             end
%         end
%         
%         % Now adjust massflow
%         if (changeMdotSign == true)
%             changeMdotSign = false;
%             d_m_dot = -1 * d_m_dot * d_m_dot_multiplier;
%             peaks = peaks + 1;
%             disp(strcat(['Crossed peak ' num2str(peaks) ', Changing direction.']));
%         end
%         
        % set new massflowrate
%         params.m_dot_IP = params.m_dot_IP + d_m_dot; 
    result.secondsToSolve = toc;
            
    end