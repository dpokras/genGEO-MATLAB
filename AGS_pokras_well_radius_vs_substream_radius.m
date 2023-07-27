clear;
close all;
clc;

date = '_22_06_23';

%% Define System Parameters
params                  = SimulationParameters; %radius, surface temp, dT/dz, etc.
params.wellCostType     = 'Baseline';
params.config  = 1;
params.S_ratio = 0; %S_ratio = m_S1/m_total
params.find_opt_S_ratio_toggle = 0;
if params.config == 1 || params.config == 2
    params.S_ratio = 1;
end
params.res_length = 7000; %m
params.n_streams = 7;

%% Solver toggle    
useMySolver = 1;

%% Collection
radii = linspace(0.2,1.0,5); %m
radius_ratios = linspace(0.25,1.0,4);
A = zeros(size(radius_ratios,2),19);
S = struct;
results = struct;
Z = zeros(size(radii,2),size(radius_ratios,2));

for i = 1:size(radii,2)
    
    %this speeds up the optimisation by bringing the starting value closer
    %to the optimal    
    filename = strcat('dpo_opt_config1_wellradius_',num2str(radii(i)),date,'m.csv');
    %% Calculate
    
    for j = 1:size(radius_ratios,2)
        
        x0 = 50;

        params.well_radius = radii(i);
        params.side_stream_radius = radii(i)*radius_ratios(j);
        
        tic;
        [result,x_opt] = AGS_solver(useMySolver,x0,params);

        if result.s_turb_in <= (1.65e3*0.99)
                disp('constraint NOT satisfied')
                con = 0;
        else
            disp('constraint satisfied')
            con = 1;
        end

        toc;                
        A(j,:) = [result.Power_net_avg/1e3,...
            result.Power_reservoir_avg/1e3,...
            result.SpCC_W_brownfield,...
            result.SpCC_Q_brownfield,...
            result.SpCC_W_greenfield,...
            result.SpCC_Q_greenfield,...
            params.well_radius,...
            radius_ratios(j),...
            params.depth,...
            params.res_length,...
            params.n_streams,...
            x_opt(1),...
            result.max_speed,...
            result.CapitalCost.C_brownfield,...
            result.CapitalCost.C_greenfield,...
            result.CapitalCost.CostSurfacePlant,...
            result.CapitalCost.CostSubsurface_brownfield,...
            result.CapitalCost.CostSubsurface_greenfield,...
            con];
        results(j).result = result;
    end

    T = array2table(A,...
        'VariableNames',{'W_net', ...
        'W_reservoir',...
        'SpCC_W_Br', ...
        'SpCC_Q_Br', ...
        'SpCC_W_Gr', ...
        'SpCC_Q_Gr', ...
        'radius',...
        'radius_ratio',...
        'depth', ...
        'res_length', ...
        'n_streams', ...
        'm_dot', ...
        'max_speed', ...
        'C_Br', ...
        'C_Gr',...
        'C_SurfacePlant',...
        'C_Subsurface_brownfield',...
        'C_Subsurface_greenfield',...
        'con satisfied'});
    writetable(T, filename);
    resultss(i).results = results;
    S(i).T = T;
end
save('nstreamVSreslengthresults.mat', 'S');