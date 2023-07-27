params = SimulationParameters;

params.time_years       = 30;
params.wellCostType     = 'Baseline';
params.depth = 6000;
params.m_dot = 40;
params.n_streams = 6;
params.side_stream_radius = 0.1;
params.well_radius = 0.15;

%%%%%
%%%%% Current depth = 6000m
%%%%% If angle is for example: 45% then the additional depth for the
%%%%% lateral wells are L_total*sin(angle). Knowing this, the max depth for
%%%%% the last lateral well must be 6000m - L_total*sin(angle).
%%%%%

P_inj_surface = 6.5e6; %Pa
T_inj_surface = 22; %C

best = 1e6;
iter = 1;
i = 1;
test_criteria = [];
angles = [];
colors = ["black","red", "blue", "green"];
figure();
hold on;
for angle = 0
    params.angle = angle;
    angles = [angles,angle];
    for n = 10:60:190
        res_length_total = 5000*6;
        params.res_length = linspace(5000,5000+n*100, 6)*res_length_total/sum(linspace(5000,5000+n*100, 6));
        params.res_length = flip(params.res_length);
%         params.res_length = ones(1,6)*mean(params.res_length);
    
        params.depth = 6000;
%         params.depth = params.depth - max(params.res_length)*sin(params.angle*pi/180);
        depths = params.depth - params.stream_spacing*((params.n_streams-1):-1:0);
        temp_at_depth = params.T_surface_rock_C + params.dT_dz * depths; %C
    
        
        result_injWell = semi_analytic_well(P_inj_surface, T_inj_surface, params.T_surface_rock_C, -1*depths, 0, params.m_dot, params.well_radius, params);  
        result_reservoir = semi_analytic_well(result_injWell.EndPressure, result_injWell.EndTemp, temp_at_depth, 0, params.res_length, params.m_dot/params.n_streams, params.side_stream_radius, params);  
        result_prodWell = semi_analytic_well(result_reservoir.EndPressure, result_reservoir.EndTemp, temp_at_depth, params.depth, 0, params.m_dot, params.well_radius, params);
        
    
        T_prod_surface = mean([result_prodWell.EndTemp]);
        P_prod_surface = result_prodWell.EndPressure(1);    
        P_surfacePlant = Surface_plant.plant_dP(P_prod_surface, T_prod_surface, params);
        result_turbine = Surface_plant.TurbWorkfunc(P_surfacePlant, T_prod_surface, params);
        W_turbine = result_turbine.W_turbine;
        [W_cond_HX, Q_cond_HX] = Surface_plant.CondWorkfunc(result_turbine.P, result_turbine.T, params);
        W_pump = 0;
        W_net = Surface_plant.NetWorkfunc(W_turbine, W_cond_HX, W_pump);
        result_capitalCost = CapitalCost_CPG(W_turbine, 0, Q_cond_HX, W_pump, 0, params);
        result_capitalCost.CostSurfacePlant.C_pump_orc = 0;
        result_capitalCost.CostSurfacePlant.C_heatExchanger = 0;
        result_capitalCost.CostSurfacePlant.C_recuperator = 0;
        result_capitalCost.CostSurfacePlant.C_pump_prod = 0;
        result_brownfield = LCOE_Simple(result_capitalCost.C_brownfield, W_net, params);
        result_greenfield = LCOE_Simple(result_capitalCost.C_greenfield, W_net, params);
        test_criteria =[test_criteria,result_greenfield.LCOE*1e3];
        if test_criteria < best
    %         best = result_brownfield.LCOE*1e3;
            best = test_criteria;
    %         start_length_best = start_length;
        end
        iter = iter + 1;
    
    
    patchline(result_reservoir.res_length(:,1)/1e3, result_reservoir.Temp(:,1), ...
        'edgecolor', colors(i), ...
        'linewidth',1.5, ...
        'edgealpha',0.5, ...
        'facecolor', 'none',...
        'DisplayName', strcat('config', num2str(i)) ...
         );


    for k=2:6
        patchline(result_reservoir.res_length(:,k)/1e3, result_reservoir.Temp(:,k), ...
            'edgecolor', colors(i), ...
            'linewidth', 1.5, ...
            'edgealpha', 0.5, ...
            'facecolor', 'none',...
            'HandleVisibility','off'...
             );
    end
    i = i + 1;
    end
end

legend('Location','northwest');
xlabel('lateral well length [km] (max depth = 6000m)');
ylabel('Temperature [C]');

%%%%%%%%

X = categorical((1:4));
figure();
for k = 1:4
    hold on;
    h = bar(X(k),test_criteria(k));
    set(h,'facecolor',colors(k));
end
xlabel('downward angle [deg] (max depth = 6000m)');
ylabel('greenfield LCOE [$/kWh]');

%%%%%%%%

%%%%%%%%

X = categorical(angles);
figure();
for k = 1:4
    hold on;
    h = bar(X(k),test_criteria(k));
    set(h,'facecolor',colors(k));
end
xlabel('downward angle [deg] (max depth = 6000m)');
ylabel('greenfield LCOE [$/kWh]');

%%%%%%%%

% params.res_length = linspace(4000,4000+25*100, 6);
% params.res_length = flip(params.res_length);
% params.res_length = ones(1,6)*mean(params.res_length);
% 
% result_reservoir = semi_analytic_well(P_inj_downhole, T_inj_downhole, temp_at_depth, 0, params.res_length, params.m_dot/params.n_streams, params.side_stream_radius, params);  
% result_prodWell = semi_analytic_well(result_reservoir.EndPressure, result_reservoir.EndTemp, temp_at_depth, params.depth, 0, params.m_dot, params.well_radius, params);
% T_prod_surface = mean([result_prodWell.EndTemp]);
% P_prod_surface = result_prodWell.EndPressure(1);    
% P_surfacePlant = Surface_plant.plant_dP(P_prod_surface, T_prod_surface, params);
% result_turbine = Surface_plant.TurbWorkfunc(P_surfacePlant, T_prod_surface, params);
% W_turbine = result_turbine.W_turbine;
% disp(W_turbine/1e3);
% [W_cond_HX, Q_cond_HX] = Surface_plant.CondWorkfunc(result_turbine.P, result_turbine.T, params);
% W_pump = 0;
% W_net = Surface_plant.NetWorkfunc(W_turbine, W_cond_HX, W_pump);
% result_capitalCost = CapitalCost_CPG(W_turbine, 0, Q_cond_HX, W_pump, 0, params);
% result_capitalCost.CostSurfacePlant.C_pump_orc = 0;
% result_capitalCost.CostSurfacePlant.C_heatExchanger = 0;
% result_capitalCost.CostSurfacePlant.C_recuperator = 0;
% result_capitalCost.CostSurfacePlant.C_pump_prod = 0;
% result_brownfield = LCOE_Simple(result_capitalCost.C_brownfield, W_net, params);
% result_greenfield = LCOE_Simple(result_capitalCost.C_greenfield, W_net, params);

colors = ['red', 'blue', 'green', 'black']

figure();
hold on;
for i=1:6
plot(result_reservoir.res_length(:,i), result_reservoir.Temp(:,i), ...
    'Color', rand(1, 3), ...
    'LineWidth',1.5, ...
    'DisplayName',strcat('stream no. ',num2str(i)) ...
    );
end
ylim([60 190])
legend('Location','southeast');
title('reservior end temperature');
text(1000, 85, strcat('Turbine Power: ', num2str(W_turbine/1e3), ' kW'))
text(1000, 80, strcat('horizontal well cost: ', num2str((result_capitalCost.C_wells_horizontal/1e6)/(W_turbine/1e3)), ' M$/kW'))
text(1000, 75, strcat('Mean Temp: ', num2str(mean(result_reservoir.Temp(end,:))), ' C'))
text(1000, 70, strcat('Total res length: ', num2str(sum(params.res_length)/1e3), ' km'))
text(1000, 65, strcat('max depth: ', num2str(sum(params.depth(end))/1e3), ' km'))

% figure();
% hold on;
% for i=1:6
% plot(result_reservoir.res_length(:,i), result_reservoir.Enthalpy(:,i), ...
%     'Color', rand(1, 3), ...
%     'LineWidth',1.5, ...
%     'DisplayName',strcat('stream no. ',num2str(i)) ...
%     );
% end
% legend('Location','southeast');
% title('reservior end enthalpy');
% text(1000, 3.56e5, strcat('Turbine Power:  ',num2str(W_turbine/1e3), ' kW'))
% text(1000, 3.49e5, strcat('LCOE (green): ', num2str(result_capitalCost.C_wells_horizontal/1e6), ' M$'))
% text(1000, 3.42e5, strcat('Mean Enthalpy: ', num2str(mean(result_reservoir.Enthalpy(end,:))/1e3), ' kJ/kg'))
% text(1000, 3.35e5, strcat('Total res length: ', num2str(sum(params.res_length)/1e3), ' km'))
% text(1000, 3.28e5, strcat('max depth: ', num2str(sum(params.depth(end))/1e3), ' km'))

