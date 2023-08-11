clear;
close all;
clc;

date = '_04_08_23';

%% Define System Parameters
params                  = SimulationParameters; %radius, surface temp, dT/dz, etc.

params.config = 1
params.S_ratio = 0; %S_ratio = m_S1/m_total
params.find_opt_S_ratio_toggle = 0;
if params.config == 1 || params.config == 2
    params.S_ratio = 1;
end

params.well_radius = 0.4667/2; %m ~18 3/8" inner diameter
params.side_stream_radius = 0.2032/2; %m ~8" inner diameter

%% Solver toggle    
useMySolver = 1;

%% Collection
j = 5;
k = 10;
length = 11000;
A = zeros(size(j:k,2),21);

for l = 1:size(length,2)
    
    %this speeds up the optimisation by bringing the starting value closer
    %to the optimal

    params.res_length = length(l);
    
    filename = strcat('Results/dpo_opt_config1_res',num2str(length(l)),date,'.csv');
    %% Calculate
    
    for i = j:k
        
        x0 = 6*(size(i,2)-1)+8*(l-1)+20;

        params.n_streams = i;
        
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
        A(i-j+1,:) = [result.W_net/1e3,...
        result.Q_net/1e3,...
        result.Power/1e3,...
        result.Power_electric_sold_avg,...
        result.Power_heat_sold_avg,...
        result.Power_total_sold_avg,...
        result.CapitalCost.C_brownfield/result.Power_electric_sold_avg,...
        result.CapitalCost.C_brownfield/result.Power_heat_sold_avg,...
        result.CapitalCost.C_brownfield/result.Power_total_sold_avg,...
        result.CapitalCost.C_greenfield/result.Power_electric_sold_avg,...
        result.CapitalCost.C_greenfield/result.Power_heat_sold_avg,...
        result.CapitalCost.C_greenfield/result.Power_total_sold_avg,...
        params.depth,...
        params.res_length,...
        params.n_streams,...
        params.S_ratio,...
        x_opt(1),...
        result.max_speed,...
        result.CapitalCost.C_brownfield,...
        result.CapitalCost.C_greenfield,...
        con];
        results(i-j+1).result = result;
    end
T = array2table(A,...
    'VariableNames',{'W_net', ...
    'Q_net',...
    'Total_power',...
    'W_net_avg', ...
    'Q_net_avg',...
    'Total_power_avg',...
    'SpCC_W_Br', ...
    'SpCC_Q_Br', ...
    'SpCC_dH_Br', ...
    'SpCC_W_Gr', ...
    'SpCC_Q_Gr', ...
    'SpCC_dH_Gr', ...
    'depth', ...
    'res_length', ...
    'n_streams', ...
    'S_ratio',...
    'm_dot', ...
    'max_speed', ...
    'C_Br', ...
    'C_Gr',...
    'con satisfied'});
    writetable(T, filename);
end