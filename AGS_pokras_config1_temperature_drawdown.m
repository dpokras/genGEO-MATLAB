clear;
close all;
clc;

date = '_04_08_23';
%% Define System Parameters
params                  = SimulationParameters; %radius, surface temp, dT/dz, etc.

params.config = 1;
params.S_ratio = 1; %S_ratio = m_S1/m_total
params.find_opt_S_ratio_toggle = 0;
if params.config == 1 || params.config == 2
    params.S_ratio = 1;
end

params.well_radius = 0.4667/2; %m ~18 3/8" inner diameter
params.side_stream_radius = 0.2032/2; %m ~8" inner diameter
params.res_length = 7000;
params.n_streams = 7;

%% Solver toggle    
useMySolver = 1;

%% Collection
A = zeros(30,8);
S = struct;
results = struct;

filename = strcat('Results/dpo_opt_config1_temperature_drawdown',num2str(date),'.csv');
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
A = [result.Power/1e3,...
    ones(1,32)*result.Power_total_sold_avg/1e3,...
    ones(1,32)*result.CapitalCost.C_brownfield/result.Power_total_sold_avg,...
    ones(1,32)*result.CapitalCost.C_greenfield/result.Power_total_sold_avg,...
    ones(1,32)*params.depth,...
    ones(1,32)*params.res_length,...
    ones(1,32)*params.n_streams,...
    ones(1,32)*result.opt_m_dot,...
    [0.1,0.5,1:1:30]]
A = reshape(A, [32, 9])

T = array2table(A,...
    'VariableNames',{'Power_per_year', ...
    'dH_sold_avg', ...
    'SpCC_dH_Br', ...
    'SpCC_dH_Gr', ...
    'depth', ...
    'res_length', ...
    'n_streams', ...
    'm_dot', ...
    'year'});
writetable(T, filename);