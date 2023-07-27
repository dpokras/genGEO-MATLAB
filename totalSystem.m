clear all;
close all;

params = SimulationParameters;
params.time_years = 1;
%params.time_years = 30;
params.dT_dz = 0.030;
%params.fluid = 'CO2';
params.fluid = 'Water';
params.system = 'Porous';
%params.system = 'Conduction1';
%params.system = 'Conduction2';
%params.system = 'Conduction4';
%params.system = 'Conduction8';
params.optimizationMode = 'MinimizeLCOE_Brownfield';
%params.optimizationMode = 'MinimizeLCOE_Greenfield';
%params.optimizationMode = 'MaximizePower';
params.orcFluid = 'R245fa';
%params.orcFluid = 'R600a';
%params.wellFieldType = 'Doublet';
params.wellFieldType = '5spot_SharedNeighbor';
%params.wellCostType = 'Ideal';
%params.hasSurfaceGatheringSystem = false;

%set conditionally
params.transmissivity = 1e-15 * 50 * 300;
params.res_length = 707;

% parameter space
depth = 1000:100:7000;
%depth = 1000:1000:7000;
%depth = 1000:1000:8000;
res_length = [707];
%res_length = 1000:1000:10000;
logTrans = 2:0.125:7;
%logTrans = 2:1:7;
%logTrans = [4.18];
transmissivity = 1e-15 * 10.^logTrans;
%transmissivity = 1e-15 * [100, 180, 320, 560, 1000, 1800, 3200, 5600, 10000, 18000, 32000, 56000, 100000];
y_var = 'transmissivity';
%y_var = 'res_length';


simTime = datestr(now,'yyyymmdd_HHMMSS');
statusFilename = strcat(['data\status_' simTime '.txt']);
countFilename = strcat(['data\count_' simTime '.txt']);
filename = strcat(['data\data_' simTime '.xlsx']);

%Initialize output file with column headers
Output = [{'time_years'} {'depth'} {'res_length'} {'dT_dz'} {'T_surface_air_C'} {'T_surface_rock_C'} {'dT_approach'} ...
                {'transmissivity'} {'wellFieldType'} {'fluid'} {'system'} {'coolingMode'} {'orcFluid'} {'optimizationMode'} {'silicaPrecip'} ...
                {'m_dot_IP'} {'Q_fluid_IP_MW'} {'W_net_IP_MW'} {'Field_IP_Multiplier'} {'SpecificCapitalCost_brownfield_kW'} {'SpecificCapitalCost_greenfield_kW'} {'LCOE_brownfield_MWh'} {'LCOE_greenfield_MWh'} ...
                {'T_prod_surface_C'} {'dP_surface_MPa'} {'SecondsToSolve'} ...
                {'C_T_G_Field'} {'C_pump_orc_Field'} {'C_coolingTowers_Field'} {'C_heatExchanger_Field'} {'C_recuperator_Field'} {'C_pump_prod_Field'} {'C_pump_inj_Field'} ...
                {'C_plant_otherEquipment_Field'} {'C_plant_installation_Field'} {'C_plant_indirectContingency_Field'} {'C_surfacePlant_Field'} ...
                {'C_gatheringSystem_Field'} {'C_wells_production_Field'} {'C_wells_injection_Field'} {'C_wellfield_Field'} {'C_exploration_Field'} {'C_stimulation_Field'} ...
                {'C_brownfield_Field'} {'C_greenfield_Field'} ];

            
            
if (strcmp(y_var, 'transmissivity'))
    y = transmissivity;
elseif (strcmp(y_var, 'res_length'))
    y = res_length;
else
    throw(MException('totalSystem:NotImplemented','Not Implemented'));
end
   
%set up simulations
for i = 1:size(depth,2)
    for j = 1:size(y, 2)
        x = (i-1)*size(y,2) +j;
        param(x) = params;
        param(x).depth = depth(i);

        if (strcmp(y_var, 'transmissivity'))
            param(x).transmissivity = transmissivity(j);
            param(x).res_length = res_length(1);
        elseif (strcmp(y_var, 'res_length'))
            param(x).transmissivity = transmissivity(1);
            param(x).res_length = res_length(j);
        else
            throw(MException('totalSystem:NotImplemented','Not Implemented'));
        end
    end
end

totalIterations = size(param,2);

% The following is for parallel computing and requires the parallel
% computing toolbox be installed. Everythign up the the parfor line may be
% commented and replaced with only a for loop to run on only a single core.
poolobj = gcp('nocreate');
if (isempty(poolobj))
    %create a parallel pool with 6 cores
    poolobj = parpool(6);
end
parfor x = 1:totalIterations
%for x = 1:totalIterations
    try
        p = param(x);
       
        % try to get iteration  number from file
        try
            fid = fopen(countFilename,'r');
            val = fscanf(fid,'%s');
            iteration = str2double(val);
            if (isnan(iteration))
                iteration = 0;
            end
            fclose(fid);
        catch
            iteration = 0;
        end
        
        % add one
        iteration = iteration + 1;
        %write back
        try
            fid = fopen(countFilename,'wt');
            fprintf(fid, '%s', num2str(iteration));
            fclose(fid);
        catch
        end
        
        %display and write status text
        statusText = strcat([
            'Iteration ' num2str(iteration) ' of ' num2str(totalIterations) ...
            'x=' num2str(x) ...
            ' (' num2str(iteration/totalIterations*100, '%.1f') '%). ' ...
            'fluid=' p.fluid ', ' ...
            'orcFluid=' p.orcFluid ', ' ...
            'depth=' num2str(p.depth) ', ' ...
            'transmissivity=' num2str(p.transmissivity) ', ' ...
            'res_length=' num2str(p.res_length) '.'
            ]);
        disp(statusText);
        %Write status to disk
        try
            fid = fopen(statusFilename,'wt');
            fprintf(fid, '%s', statusText);
            fclose(fid);
        catch
        end
    
    
    
        %run simulation
        result = total_analytic_system_optmdot(p);

        W_net_IP_MW = result.W_net_IP / 1e6;
        Q_fluid_IP_MW = result.Q_fluid_IP / 1e6;
        dP_turbine_MPa = result.dP_turbine / 1e6;
        SpecificCapitalCost_brownfield_kW = result.SpecificCapitalCost_brownfield * 1e3;
        SpecificCapitalCost_greenfield_kW = result.SpecificCapitalCost_greenfield * 1e3;
        LCOE_brownfield_MWh = result.LCOE_brownfield * 1e6;
        LCOE_greenfield_MWh = result.LCOE_greenfield * 1e6;

        %Output(x) = result.W_net_IP;
        Output(x+1,:) = [p.time_years p.depth p.res_length p.dT_dz p.T_surface_air_C p.T_surface_rock_C p.dT_approach ...
            p.transmissivity p.wellFieldType {p.fluid} {p.system} {p.coolingMode} {p.orcFluid} {p.optimizationMode} {p.silicaPrecip} ...
            result.m_dot_IP ...
            Q_fluid_IP_MW ...
            W_net_IP_MW ...
            result.Field_IP_Multiplier ...
            SpecificCapitalCost_brownfield_kW ...
            SpecificCapitalCost_greenfield_kW ...
            LCOE_brownfield_MWh ...
            LCOE_greenfield_MWh ...
            result.T_prod_surface_C ...
            dP_turbine_MPa ...
            result.secondsToSolve ...
            result.CapitalCost.CostSurfacePlant.C_T_G ...
            result.CapitalCost.CostSurfacePlant.C_pump_orc ...
            result.CapitalCost.CostSurfacePlant.C_coolingTowers ...
            result.CapitalCost.CostSurfacePlant.C_heatExchanger ...
            result.CapitalCost.CostSurfacePlant.C_recuperator ...
            result.CapitalCost.CostSurfacePlant.C_pump_prod ...
            result.CapitalCost.CostSurfacePlant.C_pump_inj ...
            result.CapitalCost.CostSurfacePlant.C_plant_otherEquipment ...
            result.CapitalCost.CostSurfacePlant.C_plant_installation ...
            result.CapitalCost.CostSurfacePlant.C_plant_indirectContingency ...
            result.CapitalCost.C_surfacePlant ...
            result.CapitalCost.C_gatheringSystem ...
            result.CapitalCost.C_wells_production ...
            result.CapitalCost.C_wells_injection ...
            result.CapitalCost.C_wellfield ...
            result.CapitalCost.C_exploration ...
            result.CapitalCost.C_stimulation ...
            result.CapitalCost.C_brownfield ...
            result.CapitalCost.C_greenfield ];
    catch ME
        disp(strcat(['x=' num2str(i) ', y=' num2str(j) ', exception: ' ME.message]));

        % Write existing data
        %xlswrite(filename, Output);

        % End
        rethrow(ME);

    end
end

% Write Results
xlswrite(filename, Output);
disp(strcat(['Saved results to ' filename]));

% Play sound for finish
load train;
%load handel;
%load gong;
sound(y,Fs);


