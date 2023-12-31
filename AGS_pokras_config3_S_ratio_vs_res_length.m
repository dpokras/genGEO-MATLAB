clear;
close all;
clc;

date = '_26_06_23';
%% Define System Parameters
params                  = SimulationParameters; %radius, surface temp, dT/dz, etc.

params.config = 3;
params.S_ratio = 0; %S_ratio = m_S1/m_total
params.find_opt_S_ratio_toggle = 0;
if params.config == 1 || params.config == 2
    params.S_ratio = 1;
end

params.well_radius = 0.4667/2; %m ~18 3/8" inner diameter
params.side_stream_radius = 0.2032/2; %m ~8" inner diameter
params.res_length = 7000; %m
params.n_streams = 7;

%% Solver toggle    
useMySolver = 1;

%% Collection
ratios = linspace(0,1.0,10);
length = linspace(5000,10000,6);
A = zeros(size(ratios,2),15);
S = struct;
results = struct;

for i = 1:size(ratios,2)        
    filename = strcat('dpo_opt_config3_S_ratio_',num2str(params.S_ratio),'_vs_res_length',num2str(date),'_.csv');
    
    for l = 1:size(length,2)
        params.res_length = length(l);
        params.S_ratio = ratios(i);
        %% Calculate
        x0 = 70;
        
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
        A(i,:) = [result.Power_total_sold_avg/1e3,...
            result.CapitalCost.C_brownfield/result.Power_total_sold_avg,...
            result.CapitalCost.C_greenfield/result.Power_total_sold_avg,...
            params.S_ratio,...
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
        results(i).result = result;
    end
    T = array2table(A,...
        'VariableNames',{'dH_sold_avg', ...
        'SpCC_dH_Br', ...
        'SpCC_dH_Gr', ...
        'S_ratio',...
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
    S(i).result = result;
end
save('S_ratios_optimised.mat', 'S');