function result = ORC_Cycle(T_in_C, params)
    
    persistent Tboil;
    persistent Pboil;
    persistent FluidType;
    persistent OptimizationMode;
    
    if (strcmp(params.orcFluid,'R245fa'))
    elseif (strcmp(params.orcFluid,'R600a'))
    else
        throw(MException('ORC_Cycle:unknownOrcFluidType','Unknown ORC Fluid Type'));
    end
    
    if (strcmp(params.optimizationMode,'MaximizePower'))
        optMode = 'maxPower';
    elseif (strcmp(params.optimizationMode,'MinimizeLCOE_Brownfield'))
        optMode = 'minCost';
    elseif (strcmp(params.optimizationMode,'MinimizeLCOE_Greenfield'))
        optMode = 'minCost';
    else
        throw(MException('ORC_Cycle:unknownOrcOptimizationMode','Unknown ORC OptimizationMode'));
    end
    
    if (strcmp(params.orcCycleType,'Subcritical'))
        if (isempty(Tboil) || ~strcmp(params.orcFluid,FluidType) || ~strcmp(params.optimizationMode,OptimizationMode))
            FluidType = params.orcFluid;
            OptimizationMode = params.optimizationMode;
            Tboil = readmatrix(fullfile('data', strcat(['ORC_Tboil_optimum_' optMode '_' params.orcFluid '.csv'])));
        end
        T_boil_C = interp1(Tboil(:,1),Tboil(:,2),T_in_C);
        dT_pinch = interp1(Tboil(:,1),Tboil(:,3),T_in_C);

        result = ORC_Cycle_Tboil(T_in_C, T_boil_C, dT_pinch, params);
        result.q_recuperator = 0;
        result.dT_LMTD_recuperator = NaN;
            
    elseif (strcmp(params.orcCycleType,'Supercritical'))
        disp('Warning: Using supercritical cycle');
        % minimum of 170C required
        if (isempty(Pboil) || ~strcmp(params.orcFluid,FluidType) || ~strcmp(params.optimizationMode,OptimizationMode))
            FluidType = params.orcFluid;
            OptimizationMode = params.optimizationMode;
            Pboil = readmatrix(fullfile('data', strcat(['ORC_Pboil_optimum_' optMode '_' params.orcFluid '.csv'])));
        end
        P_boil_Pa = interp1(Pboil(:,1),Pboil(:,2),T_in_C);

        result = ORC_Cycle_Supercrit_Pboil(T_in_C, P_boil_Pa, params);
        result.q_preheater = 0;
        result.dT_LMTD_preheater = NaN;
        
    else
        throw(MException('ORC_Cycle:unknownCycleType','Unknown Cycle Type'));
    end
end

