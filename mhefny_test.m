clear; clc
close all;

%%
params = SimulationParameters;
params.depth = 2500;
params.transmissivity = 3.0*1e-13 * 100;
params.thickness = 100;
%params.time_years = 1;
params.system = 'Porous';
params.optimizationMode = 'MinimizeLCOE_Brownfield';
params.CapacityFactor = 0.9;

params.dT_dz = 0.038;

params.orcFluid = 'R245fa';
%params.orcFluid = 'R600a';
params.wellFieldType = 'Doublet';
params.wellCostType = 'Ideal';
params.hasSurfaceGatheringSystem = false;

params.Lifetime = 10;

%% find optimal flowrate
params.time_years = 10;
params.fluid = 'CO2';
m_dot_co2 = total_analytic_system_optmdot(params).m_dot_IP;
disp(strcat(['Found CO2 mass flowrate of ' num2str(m_dot_co2)]));

%%
%params.fluid = 'Water';
%m_dot_water = total_analytic_system_optmdot(params).m_dot_IP;
%disp(strcat(['Found Water mass flowrate of ' num2str(m_dot_water)]));
%%
%find output at all years
t = 1:1:params.time_years;
for i = 1:size(t, 2)
    params.time_years = t(i);
    disp(strcat(['Trying a time of ' num2str(params.time_years)]));

    %first co2
    params.fluid = 'CO2';
    params.m_dot_IP = m_dot_co2;
    result_co2 = total_analytic_system(params);
    result_co2.time_years = params.time_years;
    result_co2.m_dot_IP = m_dot_co2;
    Output_CO2(i) = result_co2;
   disp(strcat(['CO2 Power: ' num2str(result_co2.W_net_IP/1e6) ' MWe']));

    %water
%    params.fluid = 'Water';
%    params.m_dot_IP = m_dot_water;
%    result_water = total_analytic_system(params);
%    result_water.time_years = params.time_years;
%    result_water.m_dot_IP = m_dot_water;
%    Output_Water(i) = result_water;
    
%    disp(strcat(['CO2 Power: ' num2str(result_co2.W_net_IP/1e6) ' MWe, Water Power: ' num2str(result_water.W_net_IP/1e6) ' MWe']));
end
%%
Output_CO2 = struct2table(Output_CO2);
%Output_Water = struct2table(Output_Water);
%writetable(Output_CO2,'Output_CO2.xlsx','Sheet',1);

%%

subplot(3,1,1);
%hold on
yyaxis left;
plot(Output_CO2.time_years,Output_CO2.T_prod_surface_C,'DisplayName','CO2 Prod Temp');
%plot(Output_Water.time_years,Output_Water.T_prod_surface_C,'DisplayName','Water Prod Temp');
ylabel('Production Temp [C]');
yyaxis right;
plot(Output_CO2.time_years,Output_CO2.m_dot_IP,'DisplayName','CO2 Mass Flow');
%plot(Output_Water.time_years,Output_Water.m_dot_IP,'DisplayName','Water Mass Flow');
ylabel('Opt. Mass Flowrate [kg/s]');
xlabel('Time [years]');
grid on;
legend('Location','southeast');

%%
% Mdot vs Power
subplot(3,1,2);
hold on;
plot(Output_CO2.time_years,2*Output_CO2.W_net_IP/1e6,'DisplayName','CO2 Net Power');
%plot(Output_CO2.time_years,2*Output_Water.W_net_IP/1e6,'DisplayName','Water Net Power');
hold off;
xlabel('Time [years]');
ylabel('Power [MWe]');
ylim([0 inf]);
legend('Location','northwest');
grid on;


% Time vs LCOE
subplot(3,1,3);
hold on;
plot(Output_CO2.time_years,Output_CO2.LCOE_greenfield*1e6,'DisplayName','CO2 LCOE');
%plot(Output_Water.time_years,Output_Water.LCOE_greenfield*1e6,'DisplayName','Water LCOE');
hold off;
xlabel('Time [years]');
ylabel(strcat(['LCOE [' num2str(params.costYear) '$/MW-h]']));
legend();
grid on;
%%
%set position
set(gcf,'Position',[50 50 600 900]);

fileName = strcat(['images/mhefny_manyTime_plots.png']);
saveas(gcf,fileName);
