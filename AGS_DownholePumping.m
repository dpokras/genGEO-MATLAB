function result = AGS_DownholePumping(params)

    %% Set starting parameters
    % Get condensation pressure
    T_condensation = params.T_surface_air_C + params.dT_approach;
    P_condensation = CoolProp('P', 'T', T_condensation+273.15, 'Q', 0, params.fluid) + 50e3;
    P_in = P_condensation % Pa, pump inlet pressure
    T_in = T_condensation % deg_C, pump inlet temperature
    temp_at_bottom = params.T_surface_rock_C + params.dT_dz * params.depth; %deg_C

    %% Algorithm
    % Start without pumping
    dP = 0; %pump pressure differencial
    dP_loop = 1;
    repeat = true;
    P_pump_outlet = P_in + dP;
    T_pump_out = T_in;

    % apply semi_analytic_well to Injection Well, until pressure of CO2 does
    % not go to 0
    while repeat == true
        disp(strcat(['Trying a injection pressure of ' num2str(P_pump_outlet) ]));
        P_pump_outlet = P_in + dP;
        h_in = CoolProp('HMASS', 'P', P_in, 'T', T_in+273.15, params.fluid);
        s_in = CoolProp('SMASS', 'P', P_in, 'T', T_in+273.15, params.fluid);
        % if dP_pump positive, it's a pump, if negative it's a throttle.
        if (dP > 0)
            h_pump_outletS = CoolProp('HMASS', 'P', P_pump_outlet, 'S', s_in, params.fluid);
            h_pump_outlet = h_in + (h_pump_outletS - h_in) / params.eta_cpg_pump;
        else
            h_pump_outlet = h_in;
        end
        W_pump_inj_IP = -1 * params.m_dot_IP * (h_pump_outlet - h_in); % work used by pump
        T_pump_out = CoolProp('T', 'P', P_pump_outlet, 'HMASS', h_pump_outlet, params.fluid) - 273.15;

        try
            %Injection well
            result_injWell = semi_analytic_well(P_pump_outlet, T_pump_out, params.T_surface_rock_C, -1*params.depth, 0, params.m_dot_IP, params.well_radius, params);
            repeat = false;
            %Lateral Well
            result_latWell = semi_analytic_well(result_injWell.EndPressure, result_injWell.EndTemp, temp_at_bottom, 0, params.res_length, params.m_dot_IP, params.well_radius, params);
            %Production Well
            result_prodWell = semi_analytic_well(result_latWell.EndPressure, result_latWell.EndTemp, temp_at_bottom, params.depth,0, params.m_dot_IP, params.well_radius, params);
            h_prod_surface = result_prodWell.EndEnthalpy;
        catch MException('SemiAnalytic:BelowZeroPressure','Below zero pressure!');
            dP = dP_loop * 10000 + dP; % Add 10 kPa to pump
            repeat = true;
        end
        % dP_pump can't be negative
        if (dP < 0 && dP_downhole > 0)
            dP = 0;
            disp('Warning: total_analytic_system_co2:dP_pump is zero!');
        end
    end


    % set z and T and P results for plotting
    result.z_inj = result_injWell.State(:,1) * (-1);
    result.T_inj = result_injWell.State(:,3);
    result.P_inj = result_injWell.State(:,2);
    result.P_lat = result_latWell.State(:,2);
    result.T_lat = result_latWell.State(:,3);
    result.z_lat = zeros(size(result.T_inj));
    for i = 1:size(result.z_lat)
        result.z_lat(i) = params.depth;
    end
    result.T_prod = result_prodWell.State(:,3);
    result.P_prod = result_prodWell.State(:,2);
    result.z_prod = flip(result_prodWell.State(:,1));
    result.T_start = result.T_inj(1);
    result.P_start = result.P_inj(1); % in MPa
    result.T_end = result_prodWell.EndTemp;
    result.P_end = result_prodWell.EndPressure/1e6; % in MPa
    result.dP_pump = dP;


    % set power results
    result.Q_well = result_injWell.Heat + result_latWell.Heat + result_prodWell.Heat;
    result.EndEnthalpy = h_prod_surface; % enthalpy of fluid out of production well
    result.PumpOutEnthalpy = h_pump_outlet;
    result.W_pump = W_pump_inj_IP;
    result.Q_thermal = params.m_dot_IP * (h_prod_surface - h_pump_outlet); %thermal energy collected in all wellbores


end
