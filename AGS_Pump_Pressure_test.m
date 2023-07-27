clear;
close all;
clc;
tic;

%% Define parameters for base case
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
params.res_length = 5000; % in meters
params.depth = 3500; % in meters
params.system = 'Conduction4'; % Horizontal pipe section consists of 4 laterals
params.wellFieldType = 'Doublet';    
L = 8*params.res_length + 4*params.depth;

%% Define System Parameters
m_dot_IP_max = 100; % Input maximum flowrate
params.time_years = 30;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';
%% Initialize Vectors
P_difference = zeros(1,m_dot_IP_max+1); % Initialize dP vector
m = linspace(0,m_dot_IP_max,m_dot_IP_max+1);
m = fliplr(m);
%% For Each Fluid Type and for many m-dots, find Power, Specific Power, LCOE, Relevant Pressures 
params.fluid = "Water";
for i = 1:m_dot_IP_max+1
    params.m_dot_IP = m(i)
    if params.m_dot_IP ~= 0
        result = total_analytic_system_water(params);           
        P_difference(1,i) = result.dP_pump / 1e6;
    end
end
%% Plot for Pressure for Pump (dP_Pump)
figure(1)
    plot(m,P_difference(1,:),'-b')
    ylabel('Pressure Difference [MPa]') % y-axis label
    xlim([0,m_dot_IP_max]) % Specifying limits on x-axis
    xlabel('Mass Flowrate[kg s^{-1}') % x-axis label
    title('Pressure Difference for Varying Mass Flowrates')
    legend('Pump Pressure difference Water','Turbine Pressure difference CO_2', 'Pump Pressure difference CO_2')
    grid on
    hold on
saveas(figure(1),'images\Pump_Pressure_test.fig');
toc;