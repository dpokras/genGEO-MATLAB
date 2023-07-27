% clear;
% close all;
% clc;

%% Define System Parameters
params                  = SimulationParameters; %radius, surface temp, dT/dz, etc.

params.time_years       = 30;
params.wellCostType     = 'Baseline';
% params.wellCostType     = 'Ideal';

%% Calculate
useMySolver = 1;
x0 = [6.99e3
    5e3
    11
    81
    0.16154
    0.08];
if ~useMySolver

    result = total_analytic_system_co2(wrapMyParameters(params,x0));
    disp(['The result is: [' num2str(result.W_net/1e3) '] kW']);
    disp(['LCOE Brownfield: [' num2str(result.LCOE_brownfield*1e6) '] 2019$/MWh']);
    disp(['LCOE Greenfield: [' num2str(result.LCOE_greenfield*1e6) '] 2019$/MWh']);
    
    disp(['Depth: [' num2str(x0(1)) '] m']);
    disp(['Lateral Length: [' num2str(x0(2)) '] m']);
    disp(['no of side streams: [' num2str(x0(3)) ']']);
    disp(['Mass Flow: [' num2str(x0(4)) '] kg/s']);
    disp(['Pipe Radius: [' num2str(x0(5)) '] m']);
    disp(['Side Stream Radius: [' num2str(x0(6)) '] m']);
    disp(['max speed: [' num2str(result.max_speed) ']'])
else % useMySolver

    myCostFunc = @(x)costfunc(total_analytic_system_co2(wrapMyParameters(params,x)));
%     myNonLinConstraints = @(x)contrfunc(total_analytic_system_co2(wrapMyParameters(params,x)),params);

    myOptions = optimoptions('fmincon','Display','iter',...
        'Algorithm','interior-point',...
        'MaxIterations',100);

%     x_opt = fmincon( ...
%         myCostFunc, ...
%         x0, ...
%         [],[], ...
%         [], [], ...
%        [6000;200;1;10;0.10;0.05],[8000;8000;30;inf;0.5;0.5],myNonLinConstraints,myOptions);
%     x_opt(3) = round(x_opt(3));

        x_opt = fmincon( ...
        myCostFunc, ...
        x0, ...
        [],[], ...
        [], [], ...
       [6000;200;1;10;0.10;0.05],[8000;8000;30;inf;0.5;0.5],[],myOptions);
    x_opt(3) = round(x_opt(3));

    result_opt = total_analytic_system_co2(wrapMyParameters(params,x_opt));
    disp(['Net Power: [' num2str(result_opt.W_net/1e3) '] kW']);
    disp(['LCOE Brownfield: [' num2str(result_opt.LCOE_brownfield*1e6) '] 2019$/MWh']);
    disp(['LCOE Greenfield: [' num2str(result_opt.LCOE_greenfield*1e6) '] 2019$/MWh']);
    
    disp(['Depth: [' num2str(x_opt(1)) '] m']);
    disp(['Lateral Length: [' num2str(x_opt(2)) '] m']);
    disp(['no of side streams: [' num2str(x_opt(3)) ']']);
    disp(['Mass Flow: [' num2str(x_opt(4)) '] kg/s']);
    disp(['Pipe Radius: [' num2str(x_opt(5)) '] m']);
    disp(['Side Stream Radius: [' num2str(x_opt(6)) '] m']);
    disp(['max speed: [' num2str(result_opt.max_speed) ']'])

end



