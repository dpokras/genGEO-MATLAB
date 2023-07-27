function result = ORC_Cycle_Tboil(T_in_C, T_boil_C, dT_pinch, params)

    % heat/power is output as specific heat and specific work. To find the
    % actual power, multiply by the flowrate of geofluid through the
    % system.
    
    if (isnan(T_in_C))
        throw(MException('ORC_Cycle_Tboil:T_in_NaN','ORC input temperature is NaN!'));
    end
    if (isnan(T_boil_C))
        throw(MException('ORC_Cycle_Tboil:T_boil_NaN','ORC boil temperature is NaN!'));
    end
    
    T_crit_C = CoolProp('TCRIT', "", 0, "", 0, params.orcFluid) - 273.15;
    if (T_boil_C > T_crit_C)
        throw(MException('ORC_Cycle_Tboil:Tboil_Too_Large',strcat(['Boiling temperature above critical point of ' num2str(T_crit_C)])));
    end
    
    %dT_pinch = params.dT_orc_pinch;
    if (dT_pinch <=0)
        throw(MException('ORC_Cycle_Tboil:dT_pinch_Negative',strcat(['dT_pinch is negative!'])));
    end
       
    if (T_in_C < T_boil_C + dT_pinch)
        throw(MException('ORC_Cycle_Tboil:Tboil_Too_Large',strcat(['Boiling temperature of ' num2str(T_boil_C) ' is greater than input temp of ' num2str(T_in_C) ' less pinch dT of ' num2str(params.dT_orc_pinch)])));
    end
    
    T_boil_max = MaxSubcritORCBoilTemp(params);
    %T_boil_max = 123.8;
    if (T_boil_C > T_boil_max)
        throw(MException('ORC_Cycle_Tboil:Tboil_Too_Large',strcat(['Boiling temperature of ' num2str(T_boil_C) ' is greater than maximum allowed of ' num2str(T_boil_max)])));
    end
    
    T_condense_C = params.T_surface_air_C + params.dT_approach;
    
    %State 1 (Condenser -> Pump)
    %saturated liquid
    T_C(1) = T_condense_C;
    h(1) = CoolProp('HMASS', 'T', T_C(1)+273.15, 'Q', 0, params.orcFluid);
    s(1) = CoolProp('SMASS', 'T', T_C(1)+273.15, 'Q', 0, params.orcFluid);
    P(1) = CoolProp('P', 'T', T_C(1)+273.15, 'Q', 0, params.orcFluid);
    
    %State 6 (Desuperheater -> Condenser)
    %saturated vapor
    T_C(6) = T_C(1);
    h(6) = CoolProp('HMASS', 'T', T_C(6)+273.15, 'Q', 1, params.orcFluid);
    s(6) = CoolProp('SMASS', 'T', T_C(6)+273.15, 'Q', 1, params.orcFluid);
    P(6) = P(1);
    
    %State 3 (Preheater -> Boiler)
    %saturated liquid
    T_C(3) = T_boil_C;
    h(3) = CoolProp('HMASS', 'T', T_C(3)+273.15, 'Q', 0, params.orcFluid);
    s(3) = CoolProp('SMASS', 'T', T_C(3)+273.15, 'Q', 0, params.orcFluid);
    P(3) = CoolProp('P', 'T', T_C(3)+273.15, 'Q', 0, params.orcFluid);

    %State 4 (Boiler -> Turbine)
    %saturated vapor
    T_C(4) = T_C(3);
    h(4) = CoolProp('HMASS', 'T', T_C(4)+273.15, 'Q', 1, params.orcFluid);
    s(4) = CoolProp('SMASS', 'T', T_C(4)+273.15, 'Q', 1, params.orcFluid);
    P(4) = P(3);

    %State 5 (Turbine -> Desuperheater)
    P(5) = P(1);
    h_5s = CoolProp('HMASS', 'P', P(5), 'SMASS', s(4), params.orcFluid);
    h(5) = h(4) - params.eta_orc_turbine*(h(4)-h_5s);
    T_C(5) = CoolProp('T', 'HMASS', h(5), 'P', P(5), params.orcFluid) - 273.15;
    s(5) = CoolProp('SMASS', 'HMASS', h(5), 'P', P(5), params.orcFluid);

    %State 2 (Pump -> Preheater)
    P(2) = P(3);
    h_2s = CoolProp('HMASS', 'P', P(2), 'SMASS', s(1), params.orcFluid);
    h(2) = h(1) - ((h(1)-h_2s)/params.eta_orc_pump);
    T_C(2) = CoolProp('T', 'HMASS', h(2), 'P', P(2), params.orcFluid) - 273.15;
    s(2) = CoolProp('SMASS', 'HMASS', h(2), 'P', P(2), params.orcFluid);

    %Calculate orc heat/work
    w_pump_orc = h(1)-h(2);
    q_preheater_orc = -1*(h(2)-h(3));
    q_boiler_orc = -1*(h(3)-h(4));
    w_turbine_orc = h(4)-h(5);
    q_desuperheater_orc = -1*(h(5)-h(6));
    q_condenser_orc = -1*(h(6)-h(1));
    
    dP_pump_orc = P(2)-P(1);
    P_boil = P(3);

    % Cooling Tower Parasitic load
    dT_range = T_C(5)-T_C(6);
    [f_cooling, f_condensing] = ParasiticPowerFraction_CoolingTower(params.T_surface_air_C, params.dT_approach, dT_range, params.coolingMode);
    w_cooler_orc = q_desuperheater_orc * f_cooling;
    w_condenser_orc = q_condenser_orc * f_condensing;

    %water (assume pressure 100 kPa above saturation)
    P_sat = CoolProp('P', 'T', T_in_C + 273.15, 'Q', 0, 'Water');
    cp = CoolProp('CPMASS', 'T', T_in_C + 273.15, 'P', P_sat + 100e3, 'Water');
    %Water state 11, inlet, 12, mid, 13 exit
    T_C(11) = T_in_C;
    T_C(12) = T_boil_C + dT_pinch;
    %mdot_ratio = mdot_orc / mdot_water
    mdot_ratio = cp * (T_C(11) - T_C(12)) / q_boiler_orc;
    T_C(13) = T_C(12) - mdot_ratio * q_preheater_orc / cp;

    % check that T_C(13) isn't below pinch constraint
    if (T_C(13) < T_C(2) + dT_pinch)
        %pinch constraint is here, not at 12
        % outlet is pump temp plus pinch
        T_C(13) = T_C(2) + dT_pinch;
        R = q_boiler_orc / (q_boiler_orc + q_preheater_orc);
        T_C(12) = T_C(11) - (T_C(11)-T_C(13))*R;
        mdot_ratio = cp * (T_C(11) - T_C(12)) / q_boiler_orc;
    end

    %Calculate water heat/work
    w_pump = mdot_ratio * w_pump_orc;
    q_preheater = mdot_ratio * q_preheater_orc;
    q_boiler = mdot_ratio * q_boiler_orc;
    w_turbine = mdot_ratio * w_turbine_orc;
    q_desuperheater = mdot_ratio * q_desuperheater_orc;
    q_condenser = mdot_ratio * q_condenser_orc;
    w_cooler = mdot_ratio * w_cooler_orc;
    w_condenser = mdot_ratio * w_condenser_orc;

    w_net = w_turbine + w_pump + w_cooler + w_condenser;
    w_net_orc = w_turbine_orc + w_pump_orc + w_cooler_orc + w_condenser_orc;
    
    % Calculate temperatures
    dT_range_CT = T_C(5) - T_C(6);
    dT_A_p = T_C(13) - T_C(2);
    dT_B_p = T_C(12) - T_C(3);
    if (dT_A_p == dT_B_p)
        dT_LMTD_preheater = dT_A_p;
    else
        dT_LMTD_preheater = (dT_A_p - dT_B_p) / log(dT_A_p/dT_B_p);
    end
    dT_A_b = T_C(12) - T_C(3);
    dT_B_b = T_C(11) - T_C(4);
    if (dT_A_b == dT_B_b)
        dT_LMTD_boiler = dT_A_b;
    else
        dT_LMTD_boiler = (dT_A_b - dT_B_b) / log(dT_A_b/dT_B_b);
    end
    
    T_out_C = T_C(13);
    
    
    result.q_preheater = q_preheater;
    result.q_boiler = q_boiler;
    result.w_turbine = w_turbine;
    result.q_desuperheater = q_desuperheater;
    result.q_condenser = q_condenser;
    result.w_pump = w_pump;
    result.w_cooler = w_cooler;
    result.w_condenser = w_condenser;
    result.w_net = w_net;
    result.dT_range_CT = dT_range_CT;
    result.dT_LMTD_preheater = dT_LMTD_preheater;
    result.dT_LMTD_boiler = dT_LMTD_boiler;
    result.T_out_C = T_out_C;
    result.T_condense_C = T_condense_C;
    result.T_crit_C = T_crit_C;
    result.T_boil_C = T_boil_C;
    result.dT_pinch = dT_pinch;
    result.mdot_ratio = mdot_ratio;
    result.f_cooling = f_cooling;
    result.f_condensing = f_condensing;
    result.dP_pump_orc = dP_pump_orc;
    result.P_boil = P_boil;
    
    result.T_C = T_C;
    result.P = P;
    result.h = h;
    result.s = s;
end

