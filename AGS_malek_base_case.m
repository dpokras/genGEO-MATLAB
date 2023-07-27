% clear;
% close all;clc

%% Initials
% Select Fluid Type
fluid_Id                = 2;
fluid_All               = ["Water","CO2"];

%% Define System Parameters
params                  = SimulationParameters; %radius, surface temp, dT/dz, etc.

params.wellFieldType    = 'Doublet';
params.time_years       = 30;
params.optimizationMode = 'MinimizeLCOE_Greenfield'; % MaximizePower MinimizeLCOE_Greenfield
params.orcFluid         = 'R245fa';
params.depth            = 8000;
params.dT_dz            = 0.035;

%% Base Case
params.dT_dz            = 0.035;
params.system           = "Conduction4";
params.n_streams        = 6;
% params.system = "Conduction8"; % params.res_length = 4174.1908965;

params.res_length       = 2091.5;
params.well_radius      = 0.1;
% params.well_radius = 0.244475/2; % params.res_length = 3899.643335;

% params.wellCostType     = 'Baseline';
params.wellCostType     = 'Ideal';
params.fluid            = fluid_All(fluid_Id);
params.m_dot_IP         = 43.798;

%% Calculate
useMySolver = 1;

if ~useMySolver
    result = total_analytic_system_co2(params);
    Power  = result.W_net_IP/1e3; % in kW
    LCOE   = result.LCOE_brownfield*1e6; %in 2019$/MWh
else % useMySolver
    x0 = [7.99e3
        1.8e3
        5
        40
        0.1];

    myCostFunc = @(x)costfunc(total_analytic_system_co2(wrapMyParameters(params,x)));
    myNonLinConstraints = @(x)contrfunc(total_analytic_system_co2(wrapMyParameters(params,x)),params);

    myOptions = optimoptions('fmincon','Display','iter', ...
        'MaxIterations',100);

    x_opt = fmincon( ...
        myCostFunc, ...
        x0, ...
        [],[], ...
        [], [], ...
       [1e3;500;1;10;0.05],[8000;7000;20;inf;0.25],myNonLinConstraints,myOptions);
    x_opt(3) = round(x_opt(3))
    x_opt
end

return