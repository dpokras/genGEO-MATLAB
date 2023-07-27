clear;
close all;
clc;

% date = '_21_06_23';

%% Define System Parameters
params                  = SimulationParameters; %radius, surface temp, dT/dz, etc.

params.time_years       = [1,5,15,30];
params.wellCostType     = 'Baseline';
% params.wellCostType     = 'Ideal';
params.find_opt_S_ratio_toggle = 0;
params.S_ratio = 0.1; %S_ratio = m_S1/m_total
params.config = 3;
if params.config == 1 || params.config == 2
    params.S_ratio = 1;
end

params.well_radius = 0.4667/2; %m ~18 3/8" inner diameter
params.side_stream_radius = 0.2032/2; %m ~8" inner diameter


params.coolingMode = 'Wet';

%% Solver toggle    
useMySolver = 0;
x0 = 70;

%% Calculate

params.n_streams = 9;
params.res_length = 10000;

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
%         
% A = [result.W_net/1e3,...
%     result.SpCC_W_brownfield,...
%     result.SpCC_Q_brownfield,...
%     result.SpCC_W_greenfield,...
%     result.SpCC_Q_greenfield,...
%     params.depth,...
%     params.res_length,...
%     params.n_streams,...
%     x_opt(1),...
%     result.max_speed,...
%     result.CapitalCost.C_brownfield,...
%     result.CapitalCost.C_greenfield,...
%     con];
% 
% % T = array2table(A,...
% %     'VariableNames',{'W_net', ...
% %     'SpCC_W_Br', ...
% %     'SpCC_Q_Br', ...
% %     'SpCC_W_Gr', ...
% %     'SpCC_Q_Gr', ...
% %     'depth', ...
% %     'res_length', ...
% %     'n_streams', ...
% %     'm_dot', ...
% %     'max_speed', ...
% %     'C_Br', ...
% %     'C_Gr',...
% %     'con satisfied'});
% writematrix(result.HS1,'dpo_HS1_22_06_23.csv')
% writematrix(result.hS1,'dpo_hS1_22_06_23.csv')
% writematrix(result.hS2,'dpo_hS2_22_06_23.csv')
% writematrix(result.PS2,'dpo_PS2_22_06_23.csv')
% writematrix(result.PS1,'dpo_PS1_22_06_23.csv')
% 
% figure(1)
% hold on
% plot(result.HS1/1e6, result.PS1/1e6, 'k');
% % plot(result.HS2/1e6, result.PS2/1e6, ':k');
% scatter(result.H_points/1e6, result.P_points/1e6, 'k')
% xlabel('Enthalpy Flow [MW]');
% ylabel('Pressure [MPa]');
% hold off
% 
% figure(2)
% hold on
% plot(result.hS1/1e3, result.PS1/1e6,'k');
% plot(result.hS2/1e3, result.PS2/1e6,':k');
% scatter(result.h_points/1e3, result.P_points/1e6, 'k')
% xlabel('Enthalpy [kJ/kg]');
% ylabel('Pressure [MPa]');
% hold off
% 
% D = cat(1, result.injection.stream_depth, result.reservoir.stream_depth+result.injection.stream_depth(end,:), result.production.stream_depth+result.injection.stream_depth(end,:));
% D = cat(1,D, D(end,:));
% L = cat(1, result.injection.length, result.reservoir.length+result.injection.length(end,:), result.production.length+result.reservoir.length(end,:)+result.injection.length(end,:));
% L = cat(1,L, ones(1,params.n_streams)*L(end));
% P = cat(1, result.injection.Pressure, result.reservoir.Pressure, result.production.Pressure, result.PS1(80,:));
% T_prod_surface_flashed = CoolProp.PropsSI('T', 'P', result.PS1(80,:), 'HMASS', result.production.EndEnthalpy, params.fluid) - 273.15; %C
% T_prod_surface = ones(1, params.n_streams) *  mean(T_prod_surface_flashed);
% T = cat(1, result.injection.Temp, result.reservoir.Temp, result.production.Temp, T_prod_surface);
% 
% figure(3)
% plot(L/1e3, T, 'k');
% yyaxis left
% xlabel('Length along Well [km]');
% ylabel('Temperature [C]');
% hold on
% yyaxis right
% plot(L/1e3, P/1e6,'--k');
% ylabel('Pressure [MPa]');
% hold off
% 
% figure(4)
% plot(T, D/1e3,'k');
% xlabel('Temperature [C]');
% ylabel('Depth [km]');
% 
% figure(5);
% plot(P/1e6,D/1e3,'--k')
% xlabel('Pressure [MPa]');
% ylabel('Depth [km]');
% 