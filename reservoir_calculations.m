clc;
clear;
params = SimulationParameters;

params.time_years       = 30;
params.wellCostType     = 'Baseline';
params.depth = 6000;
params.m_dot = 70;
params.n_streams = 7;
params.well_radius = 0.4667/2; %m ~18 3/8" inner diameter
params.side_stream_radius = 0.2032/2; %m ~8" inner diameter
params.config = 1;
params.S_ratio = 1;

%%%%%
%%%%% Current depth = 6000m
%%%%% If angle is for example: 45% then the additional depth for the
%%%%% lateral wells are L_total*sin(angle). Knowing this, the max depth for
%%%%% the last lateral well must be 6000m - L_total*sin(angle).
%%%%%

P_inj_surface = 6.5e6; %Pa
T_inj_surface = 22; %C

angles = [];
lineStyles = ["-", "--", ":", "-."]; % Array of different line styles
markers = ["o", "s", "d", "^", "v", ">", "<", "p", "h"]; % Array of different markers
figure();
hold on;
i = 1;
for angle = 0:30:90
    params.angle = angle;
    angles = [angles,angle];
    for n = 60
        res_length_total = 5000*6;
        params.res_length = linspace(5000,5000+n*100, 6)*res_length_total/sum(linspace(5000,5000+n*100, 6));
        params.res_length = ones(1,params.n_streams)*mean(params.res_length);

        params.depth = 6000; 
        params.depth = params.depth - params.res_length(end)*sin(angle*pi/180);
        depths = params.depth - params.stream_spacing*((params.n_streams-1):-1:0);
        temp_at_depth = params.T_surface_rock_C + params.dT_dz * depths; %C
    
        result_injWell = semi_analytic_well(P_inj_surface, T_inj_surface, params.T_surface_rock_C, -1*depths, 0, params.m_dot, params);  
        result_reservoir = semi_analytic_well(result_injWell.EndPressure, result_injWell.EndTemp, temp_at_depth, 0, params.res_length, params.m_dot/params.n_streams, params);  
        result_prodWell = semi_analytic_well(result_reservoir.EndPressure, result_reservoir.EndTemp, temp_at_depth, params.depth, 0, params.m_dot, params);

        % Overlay lines and markers
        plot(result_reservoir.res_length(:,1)/1e3, result_reservoir.Temp(:,1), ...
            'linestyle', lineStyles(i), ...
            'marker', markers(i), ...
            'MarkerIndices', 2:2:length(result_reservoir.res_length(:,1)), ...
            'color', 'black',...
            'DisplayName', strcat(num2str(angle), '°') ...
        );

        for k=2:params.n_streams
            plot(result_reservoir.res_length(:,k)/1e3, result_reservoir.Temp(:,k), ...
                'linestyle', lineStyles(i), ...
                'marker', markers(i), ...
                'MarkerIndices', 2:2:length(result_reservoir.res_length(:,1)), ...
                'color', 'black',...
                'HandleVisibility','off'...
            );
        end
    end
    i = i + 1;
end
legend('Location','northwest');
xlabel('Reservoir Length [km] (max depth = 6000m)'); % Add an x-axis label
ylabel('Fluid Temperature [C]'); % Add a y-axis label

SpCC_W_greenfield = [];
angles = [];
i = 1;
number = 1; % Initialize configuration number
markerIndex = 1; % Initialize marker index
for angle = 0
    params.angle = angle;
    angles = [angles,angle];
    for n = 0:60:190
        res_length_total = 5000*params.n_streams;
        params.res_length = linspace(5000,5000+n*100, params.n_streams)*res_length_total/sum(linspace(5000,5000+n*100, params.n_streams));
%         params.res_length = ones(1,params.n_streams)*mean(params.res_length);

        params.depth = 6000; 
        params.depth = params.depth - params.res_length(end)*sin(angle*pi/180);
        depths = params.depth - params.stream_spacing*((params.n_streams-1):-1:0);
        temp_at_depth = params.T_surface_rock_C + params.dT_dz * depths; %C
    
        result_injWell = semi_analytic_well(P_inj_surface, T_inj_surface, params.T_surface_rock_C, -1*depths, 0, params.m_dot, params);  
        result_reservoir = semi_analytic_well(result_injWell.EndPressure, result_injWell.EndTemp, temp_at_depth, 0, params.res_length, params.m_dot/params.n_streams, params);  
        result_prodWell = semi_analytic_well(result_reservoir.EndPressure, result_reservoir.EndTemp, temp_at_depth, params.depth, 0, params.m_dot, params);

        % To calculate the power, captical cost and specific cost of
        % capital
        
        P_prod_surface_flashed = ones(1, params.n_streams) .* result_prodWell.EndPressure(:,1,:);
        T_prod_surface_flashed = CoolProp.PropsSI('T', 'P', P_prod_surface_flashed, 'HMASS', result_prodWell.EndEnthalpy, params.fluid) - 273.15; %C
        T_prod_surface = mean(T_prod_surface_flashed);
        P_prod_surface = mean(P_prod_surface_flashed);

        % pressure losses between production well and surface plant
        P_surfacePlant = Surface_plant.plant_dP(P_prod_surface, T_prod_surface, params);

        %Calculate turbine work
        if params.config == 4 || params.S_ratio == 0
            result_turbine.W_turbine = 0;
            W_turbine = 0;
            result_turbine.s_turb_in = 1.7e3;
            result_turbine.T = T_prod_surface;
            result_turbine.P = P_prod_surface;
        else
            result_turbine = Surface_plant.TurbWorkfunc(P_surfacePlant, T_prod_surface, params);
            W_turbine = result_turbine.W_turbine;
            s_turb_in = result_turbine.s_turb_in;
        end

        %Calculate the main heat exchange network
        result_HEX = Surface_plant.HEXfunc(result_turbine.P, result_turbine.T, T_prod_surface, P_prod_surface, params);
        params.S_ratio = result_HEX.S_ratio;
        Q_net = result_HEX.Q_net;

        % Calculate cooling water work
        dP_cooling_water = 2e5; %Pa
        result_cooling_pump = Surface_plant.PumpWorkfunc(dP_cooling_water, result_HEX.m_dot_water_CoolLoop, params.P_cooling_water, params.T_cooling_water, 'water', params); %W

        % Calculate Storage Compressor Work
        result_comp_unit = Surface_plant.Compfunc(params);
        W_storage_comp = result_comp_unit.W_compressor;

        % Calculate Start-Up and flow control pump
        dP_SUPump_design = 5e5;
        result_SU_pump = Surface_plant.PumpWorkfunc(dP_SUPump_design, params.m_dot, result_turbine.P, result_turbine.T, 'CO2', params); %W

        % Calculate filling pump
        % head  = 10 meters to pump fluid up to specified height. Specific
        % gravity of CO2 is 1.5189
        % reasonable m_dot is 20% of the total flowrate of the process
        dP_filling_pump = 0.0981 * (params.m_dot *0.2) * 1.5189 *1e5;
        result_filling_pump = Surface_plant.PumpWorkfunc(dP_filling_pump, (params.m_dot *0.2), 6.8e6, 15, 'CO2', params); %W

        % Calculate work for the wet cooling tower
        result_cooling_tower = result_HEX.result_cooling_tower;
        W_cooling_tower = result_cooling_tower.W_cooling_tower;

        result_pumps.result_cooling_pump = result_cooling_pump;
        result_pumps.result_SU_pump = result_SU_pump;
        result_pumps.result_filling_pump = result_filling_pump;

        % Calculate net work 
        W_net = Surface_plant.NetWorkfunc(W_turbine, W_cooling_tower, W_storage_comp, result_pumps, params);

        % Calculate the total storage capacity required for contingency;
        result_tank = Surface_plant.tank(result_injWell, result_reservoir, result_prodWell, params);
        
        %Calculate Capital Cost
        result_capitalCost = CapitalCost_CPG(result_turbine, result_HEX, result_pumps, result_cooling_tower, result_tank, result_comp_unit, params);
    
        %Calculate LCOE ($/Weh)
        %Calculate Specific Capital Cost ($/We)
        result_brownfield = LCOE_Simple(result_capitalCost.C_brownfield, W_net, Q_net, params);
        result_greenfield = LCOE_Simple(result_capitalCost.C_greenfield, W_net, Q_net, params);
    
        % costs
        result.SpCC_W_brownfield = result_brownfield.SpCC_W;
        result.SpCC_W_greenfield = result_greenfield.SpCC_W;
        SpCC_W_greenfield = [SpCC_W_greenfield, result.SpCC_W_greenfield];

        number = number + 1; % Increment configuration number
    end
    i = i + 1;
end
% Define hatch angles
hatch_angles = [15, 60, 90, 120]; 
hatch_density = [15,10,5,5];

% Create figure
figure;
hold on;
barCenters = linspace(1, length(SpCC_W_greenfield), length(SpCC_W_greenfield));
barWidth = 0.8; % This can be adjusted

for i = 1:length(SpCC_W_greenfield)
    % Create patch object for each bar
    p = patch([barCenters(i)-barWidth/2, barCenters(i)+barWidth/2, barCenters(i)+barWidth/2, barCenters(i)-barWidth/2], ...
              [2, 0, SpCC_W_greenfield(i), SpCC_W_greenfield(i)], 'w');
    % Apply hatch pattern
    hatchfill(p, 'single', hatch_angles(i), hatch_density(i)); 
    
    % Annotate the bar with its configuration number
    text(barCenters(i), SpCC_W_greenfield(i)+20, ['Configuration: ' num2str(i)], 'HorizontalAlignment', 'center');
end

% Label axes
xlabel('Configuration');
ylabel('SpCC - Greenfield [$2022/W]');

%no ticks
xticks([]);

% Set y-axis limit
ylim([0, 800]);

hold off;

SpCC_W_greenfield = [];
angles = [];
i = 1;
number = 1; % Initialize configuration number
markerIndex = 1; % Initialize marker index
for angle = 0
    params.angle = angle;
    angles = [angles,angle];
    for n = 0:60:190
        res_length_total = 5000*params.n_streams;
        params.res_length =flip(linspace(5000,5000+n*100, params.n_streams)*res_length_total/sum(linspace(5000,5000+n*100, params.n_streams)));

        params.depth = 6000; 
        params.depth = params.depth - params.res_length(end)*sin(angle*pi/180);
        depths = params.depth - params.stream_spacing*((params.n_streams-1):-1:0);
        temp_at_depth = params.T_surface_rock_C + params.dT_dz * depths; %C
    
        result_injWell = semi_analytic_well(P_inj_surface, T_inj_surface, params.T_surface_rock_C, -1*depths, 0, params.m_dot, params);  
        result_reservoir = semi_analytic_well(result_injWell.EndPressure, result_injWell.EndTemp, temp_at_depth, 0, params.res_length, params.m_dot/params.n_streams, params);  
        result_prodWell = semi_analytic_well(result_reservoir.EndPressure, result_reservoir.EndTemp, temp_at_depth, params.depth, 0, params.m_dot, params);

        % To calculate the power, captical cost and specific cost of
        % capital
        
        P_prod_surface_flashed = ones(1, params.n_streams) .* result_prodWell.EndPressure(:,1,:);
        T_prod_surface_flashed = CoolProp.PropsSI('T', 'P', P_prod_surface_flashed, 'HMASS', result_prodWell.EndEnthalpy, params.fluid) - 273.15; %C
        T_prod_surface = mean(T_prod_surface_flashed);
        P_prod_surface = mean(P_prod_surface_flashed);

        % pressure losses between production well and surface plant
        P_surfacePlant = Surface_plant.plant_dP(P_prod_surface, T_prod_surface, params);

        %Calculate turbine work
        if params.config == 4 || params.S_ratio == 0
            result_turbine.W_turbine = 0;
            W_turbine = 0;
            result_turbine.s_turb_in = 1.7e3;
            result_turbine.T = T_prod_surface;
            result_turbine.P = P_prod_surface;
        else
            result_turbine = Surface_plant.TurbWorkfunc(P_surfacePlant, T_prod_surface, params);
            W_turbine = result_turbine.W_turbine;
            s_turb_in = result_turbine.s_turb_in;
        end

        %Calculate the main heat exchange network
        result_HEX = Surface_plant.HEXfunc(result_turbine.P, result_turbine.T, T_prod_surface, P_prod_surface, params);
        params.S_ratio = result_HEX.S_ratio;
        Q_net = result_HEX.Q_net;

        % Calculate cooling water work
        dP_cooling_water = 2e5; %Pa
        result_cooling_pump = Surface_plant.PumpWorkfunc(dP_cooling_water, result_HEX.m_dot_water_CoolLoop, params.P_cooling_water, params.T_cooling_water, 'water', params); %W

        % Calculate Storage Compressor Work
        result_comp_unit = Surface_plant.Compfunc(params);
        W_storage_comp = result_comp_unit.W_compressor;

        % Calculate Start-Up and flow control pump
        dP_SUPump_design = 5e5;
        result_SU_pump = Surface_plant.PumpWorkfunc(dP_SUPump_design, params.m_dot, result_turbine.P, result_turbine.T, 'CO2', params); %W

        % Calculate filling pump
        % head  = 10 meters to pump fluid up to specified height. Specific
        % gravity of CO2 is 1.5189
        % reasonable m_dot is 20% of the total flowrate of the process
        dP_filling_pump = 0.0981 * (params.m_dot *0.2) * 1.5189 *1e5;
        result_filling_pump = Surface_plant.PumpWorkfunc(dP_filling_pump, (params.m_dot *0.2), 6.8e6, 15, 'CO2', params); %W

        % Calculate work for the wet cooling tower
        result_cooling_tower = result_HEX.result_cooling_tower;
        W_cooling_tower = result_cooling_tower.W_cooling_tower;

        result_pumps.result_cooling_pump = result_cooling_pump;
        result_pumps.result_SU_pump = result_SU_pump;
        result_pumps.result_filling_pump = result_filling_pump;

        % Calculate net work 
        W_net = Surface_plant.NetWorkfunc(W_turbine, W_cooling_tower, W_storage_comp, result_pumps, params);

        % Calculate the total storage capacity required for contingency;
        result_tank = Surface_plant.tank(result_injWell, result_reservoir, result_prodWell, params);
        
        %Calculate Capital Cost
        result_capitalCost = CapitalCost_CPG(result_turbine, result_HEX, result_pumps, result_cooling_tower, result_tank, result_comp_unit, params);
    
        %Calculate LCOE ($/Weh)
        %Calculate Specific Capital Cost ($/We)
        result_brownfield = LCOE_Simple(result_capitalCost.C_brownfield, W_net, Q_net, params);
        result_greenfield = LCOE_Simple(result_capitalCost.C_greenfield, W_net, Q_net, params);
    
        % costs
        result.SpCC_W_brownfield = result_brownfield.SpCC_W;
        result.SpCC_W_greenfield = result_greenfield.SpCC_W;
        SpCC_W_greenfield = [SpCC_W_greenfield, result.SpCC_W_greenfield];

        number = number + 1; % Increment configuration number
    end
    i = i + 1;
end
% Define hatch angles
hatch_angles = [15, -45, 0, -60]; 
hatch_density = [15,10,10,5];
configs = [1,5,6,7];

% Create figure
figure;
hold on;
barCenters = linspace(1, length(SpCC_W_greenfield), length(SpCC_W_greenfield));
barWidth = 0.8; % This can be adjusted

for i = 1:length(SpCC_W_greenfield)
    % Create patch object for each bar
    p = patch([barCenters(i)-barWidth/2, barCenters(i)+barWidth/2, barCenters(i)+barWidth/2, barCenters(i)-barWidth/2], ...
              [2, 0, SpCC_W_greenfield(i), SpCC_W_greenfield(i)], 'w');
    % Apply hatch pattern
    hatchfill(p, 'single', hatch_angles(i), hatch_density(i)); 
    
    % Annotate the bar with its configuration number
    text(barCenters(i), SpCC_W_greenfield(i)+20, ['Configuration: ' num2str(configs(i))], 'HorizontalAlignment', 'center');
end

% Label axes
xlabel('Configuration');
ylabel('SpCC - Greenfield [$2022/W]');

%no ticks
xticks([]);

% Set y-axis limit
ylim([0, 800]);

hold off;
