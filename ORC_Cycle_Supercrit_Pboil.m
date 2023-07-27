function [q_boiler, w_turbine, q_recuperator, q_desuperheater, q_condenser, w_pump, w_cooler, w_condenser, w_net, ...
        dT_range_CT, dT_LMTD_boiler, dT_LMTD_recuperator, T_out_C] ...
        = ORC_Cycle_Supercrit_Pboil(T_in_C, P_boil_Pa, T_ambient_C, dT_approach, dT_pinch, eta_pump, eta_turbine, coolingMode)

    % heat/power is output as specific heat and specific work. To find the
    % actual power, multiply by the flowrate of geofluid through the
    % system.
    
    %P_boil_Pa = 4e6;
    orcFluid = 'R245fa';
    %eta_pump = 0.9;
    %eta_turbine = 0.8;
    %dT_pinch = 5;
    
    % Critical point of R245fa
    P_crit = CoolProp('PCRIT', "", 0, "", 0, orcFluid);
    % if Pboil is below critical, throw error
    if (P_boil_Pa < P_crit)
        throw(MException('ORC_Cycle_Supercrit_Pboil:lowBoilingPressure','Boiling Pressure Below Critical Pressure'));
    end
    %the line of minimum entropy to keep the fluid vapor in turbine is
    %entropy at saturated vapor at 125C. So inlet temp must provide this
    %minimum entropy.
    s_min = CoolProp('SMASS', 'T', 125 + 273.15, 'Q', 1, orcFluid);
    T_min = CoolProp('T', 'P', P_boil_Pa, 'SMASS', s_min, orcFluid) - 273.15;
    if (T_in_C - dT_pinch < T_min)
        throw(MException('ORC_Cycle_Supercrit_Pboil:lowInletTemp',strcat(['Inlet Temp below ' num2str(T_min+dT_pinch, '%.1f') ' C for Supercritical Fluid'])));
    end
    
    T_condense_C = T_ambient_C + dT_approach;
    
    %State 1 (Condenser -> Pump)
    %saturated liquid
    T_C(1) = T_condense_C;
    h(1) = CoolProp('HMASS', 'T', T_C(1)+273.15, 'Q', 0, orcFluid);
    s(1) = CoolProp('SMASS', 'T', T_C(1)+273.15, 'Q', 0, orcFluid);
    P(1) = CoolProp('P', 'T', T_C(1)+273.15, 'Q', 0, orcFluid);
    
    %State 7 (Desuperheater -> Condenser)
    %saturated vapor
    T_C(7) = T_C(1);
    h(7) = CoolProp('HMASS', 'T', T_C(7)+273.15, 'Q', 1, orcFluid);
    s(7) = CoolProp('SMASS', 'T', T_C(7)+273.15, 'Q', 1, orcFluid);
    P(7) = P(1);

    %State 2 (Pump -> Recuperator)
    P(2) = P_boil_Pa;
    h_2s = CoolProp('HMASS', 'P', P(2), 'SMASS', s(1), orcFluid);
    h(2) = h(1) - ((h(1)-h_2s)/eta_pump);
    T_C(2) = CoolProp('T', 'HMASS', h(2), 'P', P(2), orcFluid) - 273.15;
    s(2) = CoolProp('SMASS', 'HMASS', h(2), 'P', P(2), orcFluid);
    
    %water (assume pressure 100 kPa above saturation)
    P_water = CoolProp('P', 'T', T_in_C + 273.15, 'Q', 0, 'Water') + 100e3;
    
    %Guess orc_in fluid is T_C(2)
    T_C(3) = T_C(2);
    %Water in temp is T_in_C
    T_C(11) = T_in_C;
    P(4) = P(2);
    
    dT=1;
    while (abs(dT) >= 1)
        % State 4 (Boiler -> Turbine)
        % T_C(4) = T_in_C - dT_pinch;
        % Input orc/geo heat exchanger
        [q_exchanged_orc, q_exchanged_geo, dT_LMTD_boiler, T_C(4), T_C(12), T_orc, T_geo, Q, m_dot_orc, m_dot_geo, mdot_ratio] ...
            = HeatExchanger_OptMdot(T_C(3), P(4), orcFluid, T_C(11), P_water, 'Water', dT_pinch, T_min);

        h(4) = CoolProp('HMASS', 'T', T_C(4)+273.15, 'P', P(4), orcFluid);
        s(4) = CoolProp('SMASS', 'T', T_C(4)+273.15, 'P', P(4), orcFluid);

        %State 5 (Turbine -> Recuperator)
        P(5) = P(1);
        h_5s = CoolProp('HMASS', 'P', P(5), 'SMASS', s(4), orcFluid);
        h(5) = h(4) - eta_turbine*(h(4)-h_5s);
        T_C(5) = CoolProp('T', 'HMASS', h(5), 'P', P(5), orcFluid) - 273.15;
        s(5) = CoolProp('SMASS', 'HMASS', h(5), 'P', P(5), orcFluid);

        %State 3 (Recuperator -> Boiler)
        %State 6 (Recuperator -> Desuperheater)
        P(3) = P(2);
        P(6) = P(1);
        % Assume m_dot for each fluid is 1, then output is specific heat
        % exchange
        [q_recuperator_orc, dT_LMTD_recuperator, T_C_new_3, T_C(6), T_23_recuperator, T_56_recuperator, q_recuperator_vals] ...
            = HeatExchanger(T_C(2), P(3), 1, orcFluid, T_C(5), P(6), 1, orcFluid, dT_pinch);
        h(3) = CoolProp('HMASS', 'T', T_C(3)+273.15, 'P', P(3), orcFluid);
        h(6) = CoolProp('HMASS', 'T', T_C(6)+273.15, 'P', P(6), orcFluid);
        s(3) = CoolProp('SMASS', 'T', T_C(3)+273.15, 'P', P(3), orcFluid);
        s(6) = CoolProp('SMASS', 'T', T_C(6)+273.15, 'P', P(6), orcFluid);

        dT = T_C(3) - T_C_new_3;
        T_C(3) = T_C_new_3;
    end
    
    %Calculate orc heat/work
    w_pump_orc = h(1)-h(2);
    q_boiler_orc = -1*(h(3)-h(4));
    w_turbine_orc = h(4)-h(5);
    %q_recuperator_orc = q_recuperator_orc;
    q_desuperheater_orc = -1*(h(6)-h(7));
    q_condenser_orc = -1*(h(7)-h(1));
    
    % Cooling Tower Parasitic load
    dT_range_CT = T_C(6)-T_C(7);
    [f_cooling, f_condensing] = ParasiticPowerFraction_CoolingTower(T_ambient_C, dT_approach, dT_range_CT, coolingMode);
    w_cooler_orc = q_desuperheater_orc * f_cooling;
    w_condenser_orc = q_condenser_orc * f_condensing;

    %water (assume pressure 100 kPa above saturation)
%     P_sat = CoolProp('P', 'T', T_in_C + 273.15, 'Q', 0, 'Water');
%     cp = CoolProp('CPMASS', 'T', T_in_C + 273.15, 'P', P_sat + 100e3, 'Water');
%     %Water state 11 inlet, 12 exit
%     T_C(11) = T_in_C;
%     % outlet is dT_pinch higher than the recuperator outlet
%     T_C(12) = T_C(3) + dT_pinch;
%     mdot_ratio = cp * (T_C(11) - T_C(12)) / q_boiler_orc;

    %Calculate water heat/work
    w_pump = mdot_ratio * w_pump_orc;
    q_boiler = mdot_ratio * q_boiler_orc;
    w_turbine = mdot_ratio * w_turbine_orc;
    q_recuperator = mdot_ratio * q_recuperator_orc;
    q_desuperheater = mdot_ratio * q_desuperheater_orc;
    q_condenser = mdot_ratio * q_condenser_orc;
    w_cooler = mdot_ratio * w_cooler_orc;
    w_condenser = mdot_ratio * w_condenser_orc;

    w_net = w_turbine + w_pump + w_cooler + w_condenser;
    
    % Calculate temperatures
%     dT_A_b = T_C(12) - T_C(3);
%     dT_B_b = T_C(11) - T_C(4);
%     if (dT_A_b == dT_B_b)
%         dT_LMTD_boiler = dT_A_b;
%     else
%         dT_LMTD_boiler = (dT_A_b - dT_B_b) / log(dT_A_b/dT_B_b);
%     end
    
    T_out_C = T_C(12);
end

