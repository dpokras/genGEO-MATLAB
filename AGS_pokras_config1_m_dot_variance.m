clear;
close all;
clc;

date = '_31_07_23';
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
useMySolver = 0;

%% Collection
filename = strcat('Results/dpo_opt_config1_m_dot_variance',num2str(date),'.csv');
    %% Calculate
   % Define x0 as a range from 0 to 200 with step 10
x0 = 10:10:190;

% Initialize an array to store the results
results = [];

% Loop over each value in x0
for i = 1:length(x0)
    tic;
    % Pass the current value to the solver
    [result, x_opt] = AGS_solver(useMySolver, x0(i), params);

    % Extract the required values
    power_total_sold_avg = result.Power_total_sold_avg;
    capital_cost_brownfield = result.CapitalCost.C_brownfield / power_total_sold_avg;
    capital_cost_greenfield = result.CapitalCost.C_greenfield / power_total_sold_avg;

    % Append the result to the array
    results = [results; power_total_sold_avg, capital_cost_brownfield, capital_cost_greenfield, x_opt(1)];
    toc;
end

% Convert the results to a table
resultsTable = array2table(results, 'VariableNames', {'PowerTotalSoldAvg', 'CapitalCostBrownfield', 'CapitalCostGreenfield', 'm_dot'});

% Write the table to a CSV file
writetable(resultsTable, filename);