function result = total_analytic_system_co2(params)

    if (strcmp(params.fluid,'CO2') == false)
        throw(MException('total_analytic_system_CO2:Wrong_Fluid','Wrong Fluid--Not CO2'));
    end

    % Get condensation pressure
    T_inj_surface = 25; %C
    P_inj_surface = 65e5;
    
    P_reservoir = params.depth * 1000 * 9.81;
    P_reservoir_max = params.depth * 2500 * 9.81;
    
    % If a porous system, loop until the downhole pressure is near
    % hydrostatic
    dP_threshold = 2e4; %20 kPa
    iter = 1;

    while true

%         disp(strcat(['Injection P iteration: ' num2str(iter)]));
        % enthalpy state 1 (injection)
        H1 = ones(1, params.n_streams) * CoolProp.PropsSI('HMASS', 'T', T_inj_surface+273.15, 'P',P_inj_surface, params.fluid) * params.m_dot;
        h1 = H1 / params.m_dot;
        P1 = ones(1, params.n_streams) * P_inj_surface;

        % injection well
        depths = params.depth - params.thickness*((params.n_streams-1):-1:0); %m
        temp_at_depth = params.T_surface_rock_C + params.dT_dz * depths; %C

%         disp('simulating injection well')
        result_injWell = semi_analytic_well(P_inj_surface, T_inj_surface, params.T_surface_rock_C, -1*depths, 0, params.m_dot, params);
        % enthalpy state 2 (injection well)
        H2 = result_injWell.Enthalpy * params.m_dot;
        h2 = result_injWell.Enthalpy;
        P2 = result_injWell.Pressure;

        % reservoir
        result_reservoir = semi_analytic_well(result_injWell.EndPressure, result_injWell.EndTemp, temp_at_depth, 0, params.res_length, params.m_dot/params.n_streams, params);
        % enthalpy state 3 (reservoir)
        H3 = result_reservoir.Enthalpy * params.m_dot;
        h3 = result_reservoir.Enthalpy;
        P3 = result_reservoir.Pressure;

        % if reservoir pressure greater than allowable, throw error.
        if (result_reservoir.EndPressure(end) >= P_reservoir_max)
            throw(MException('TotalAnalyticSystemCO2:ExceedsMaxReservoirPressure',strcat(['Exceeds Max Reservoir Pressure of ' num2str(P_reservoir_max/1e6,3) ' MPa!'])));
        end
    
        %production well and surface result
%         disp('simulating production well')
        result_prodWell = semi_analytic_well(result_reservoir.EndPressure, result_reservoir.EndTemp, temp_at_depth, depths, 0, params.m_dot, params);
        % enthalpy state 4 (production well)
        H4 = result_prodWell.Enthalpy * params.m_dot;
        h4 = result_prodWell.Enthalpy;
        P4 = result_prodWell.Pressure;

        %%%% Take the average temperature arriving at the production well
            %%% *** only if the stream flowrates are equal***
            %%% *** otherwise average proportionally ***
        
        %%%% Take the lowest pressure of the separate side streams and
        %%%% flash the rest
        %%%%%%%%%%
        P_prod_surface_flashed = ones(1, params.n_streams) .* result_prodWell.EndPressure(:,1,:);
        T_prod_surface_flashed = CoolProp.PropsSI('T', 'P', P_prod_surface_flashed, 'HMASS', result_prodWell.EndEnthalpy, params.fluid) - 273.15; %C
        T_prod_surface = mean(T_prod_surface_flashed);
        P_prod_surface = mean(P_prod_surface_flashed);

        % enthalpy state 4 (production end flashed)
        if isinf(any(T_prod_surface)) || isinf(any(P_prod_surface))||isnan(any(T_prod_surface)) || isnan(any(P_prod_surface))
            disp('something is wrong')
        end
        H4_F = ones(1, params.n_streams) .* CoolProp.PropsSI('HMASS', 'T', T_prod_surface+273.15, 'P',P_prod_surface, params.fluid) .* params.m_dot;
        h4_F = H4_F / params.m_dot;
        P4_F = ones(1, params.n_streams) .* P_prod_surface;

        % enthalpy state 5 (stream split; S1/2)
        H5S1 = ones(1, params.n_streams) .* CoolProp.PropsSI('HMASS', 'T', T_prod_surface+273.15, 'P',P_prod_surface, params.fluid) * params.m_dot * params.S_ratio;
        h5S1 = H5S1 / params.m_dot / params.S_ratio;
        H5S2 = ones(1, params.n_streams) .* CoolProp.PropsSI('HMASS', 'T', T_prod_surface+273.15, 'P',P_prod_surface, params.fluid) * params.m_dot * (1 - params.S_ratio);       
        h5S2 = H5S2 / (params.m_dot * (1-params.S_ratio));
        P5S1 = ones(1, params.n_streams) .* P_prod_surface;
        P5S2 = ones(1, params.n_streams) .* P_prod_surface;

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

        % enthalpy state  6(turbine outlet)
        H6 = ones(1, params.n_streams) .* CoolProp.PropsSI('HMASS', 'T', result_turbine.T+273.15, 'P',result_turbine.P, params.fluid) * params.m_dot * params.S_ratio;
        h6 = H6 / params.m_dot / params.S_ratio;
        P6 = ones(1, params.n_streams) .* result_turbine.P;

        %Calculate the main heat exchange network
        result_HEX = Surface_plant.HEXfunc(result_turbine.P, result_turbine.T, T_prod_surface, P_prod_surface, params);
        params.S_ratio = result_HEX.S_ratio;
        Q_net = result_HEX.Q_net;

        % enthalpy state  7(S1 HEX1 outlet)
        H7 = ones(1, params.n_streams) .* result_HEX.h_hex1_hot_out .* params.m_dot .* params.S_ratio;
        h7 = ones(1, params.n_streams) .* result_HEX.h_hex1_hot_out;
        P7 = ones(1, params.n_streams) .* result_HEX.P_CO2_out;
        H8 = ones(1, params.n_streams) .* result_HEX.h_hex2_hot_out .* params.m_dot .* (1-params.S_ratio);
        h8 = ones(1, params.n_streams) .* result_HEX.h_hex2_hot_out;
        P8 = ones(1, params.n_streams) .* result_HEX.P_tubeside_2;
        H9 = H8; %pressure is different P9 = P7
        h9 = h8;
        P9 = P7;
        H10 = ones(1, params.n_streams) .* CoolProp.PropsSI('HMASS', 'T', result_HEX.T_CO2_out+273.15-0.001, 'P',result_HEX.P_CO2_out, params.fluid) * params.m_dot;
        h10 = H10 / params.m_dot;
        P10 = P7;

        % Concat all Hs and Ps
        HS1 = cat(1, H1, H2, H3, H4, H4_F, H5S1, H6, H7, H10);
        hS1 = cat(1, h1, h2, h3, h4, h4_F, h5S1, h6, h7, h10);
        PS1 = cat(1, P1, P2, P3, P4, P4_F, P5S1, P6, P7, P10);        
        HS2 = cat(1, H1, H2, H3, H4, H4_F, H5S2, H8, H9, H10);
        hS2 = cat(1, h1, h2, h3, h4, h4_F, h5S2, h8, h9, h10);
        PS2 = cat(1, P1, P2, P3, P4, P4_F, P5S2, P8, P9, P10);

        H_points = cat(1,H1, H2(end,:), H3(end,:), H4_F, H5S1, H5S2, H6, H7, H8, H9, H10);
        h_points = cat(1,h1, h2(end,:), h3(end,:), h4_F, h5S1, h5S2, h6, h7, h8, h9, h10);
        P_points = cat(1,P1, P2(end,:), P3(end,:), P4_F, P5S1, P5S2, P6, P7, P8, P9, P10);

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
        
        % Changing P_condenstion to be an iteration of P until liquid condensation.

        dP_inj = result_HEX.P_CO2_out - P_inj_surface;
        P_inj_surface = result_HEX.P_CO2_out;
        T_inj_surface = result_HEX.T_CO2_out - 0.001;
        
        
        if abs(dP_inj) <= dP_threshold && iter > 1
            break
        end
        iter = iter + 1;
    end
    
%     disp(strcat(['Turbine Inlet T: ' num2str(T_prod_surface) ' C'] ));
    delta_H = H4_F - H1; %W
    delta_h = h4_F - h1; %J/kg
    
    %Calculate Capital Cost
    result_capitalCost = CapitalCost_CPG(result_turbine, result_HEX, result_pumps, result_cooling_tower, result_tank, result_comp_unit, params);

    %Calculate LCOE ($/Weh)
    %Calculate Specific Capital Cost ($/We)
    result_brownfield = LCOE_Simple(result_capitalCost.C_brownfield, W_net, Q_net, params);
    result_greenfield = LCOE_Simple(result_capitalCost.C_greenfield, W_net, Q_net, params);
    
    % set result
    % set result_wells
    result.injection = result_injWell;
    result.reservoir = result_reservoir;
    result.production = result_prodWell;
    
    % set result_plant
    result.result_turbine = result_turbine;
    result.result_HEX = result_HEX;
    result.result_cooling_pump = result_cooling_pump;
    result.result_cooling_tower = result_cooling_tower;
    result.result_tank = result_tank;
    result.W_turbine = result_turbine.W_turbine;
    result.W_net = W_net;
    result.Q_net = Q_net;
    result.dH_sold = W_net + Q_net;
    result.s_turb_in = result_turbine.s_turb_in;
    
    % set result_fluid
    result.temp_at_depth = temp_at_depth;
    result.P_reservoir = P_reservoir;
    result.P_reservoir_max = P_reservoir_max;    
    result.T_prod_surface = T_prod_surface;
    result.P_prod_surface = P_prod_surface;
    result.delta_H = delta_H;
    result.delta_h = delta_h;
    result.max_speed = params.m_dot/(pi*(params.well_radius)^2)/result_prodWell.EndDensity(1);
    
    % costs
    result.CapitalCost = result_capitalCost;
    result.SpCC_W_brownfield = result_brownfield.SpCC_W;
    result.SpCC_W_greenfield = result_greenfield.SpCC_W;
    result.SpCC_Q_brownfield = result_brownfield.SpCC_Q;
    result.SpCC_Q_greenfield = result_greenfield.SpCC_Q;
    result.SpCC_dH_brownfield = result_brownfield.SpCC_dH;
    result.SpCC_dH_greenfield = result_greenfield.SpCC_dH;
    result.LCOE_brownfield = result_brownfield.LCOE;
    result.LCOE_greenfield = result_greenfield.LCOE;

    %PH diagram
    result.HS1 = HS1;
    result.HS2 = HS2;
    result.hS1 = hS1;
    result.hS2 = hS2;
    result.PS1 = PS1;
    result.PS2 = PS2;
    result.H_points = H_points;
    result.h_points = h_points;
    result.P_points = P_points;

    result.S_ratio = result_HEX.S_ratio;
end