function [Q_exchanged, dT_LMTD, T_1_out, T_2_out, T_1, T_2, Q] = HeatExchanger(T_1_in, P_1, m_dot_1, fluid_1, T_2_in, P_2, m_dot_2, fluid_2, dT_pinch)

% Only tested where Temp1 < Temp2
if (T_1_in > T_2_in)
    %throw(MException('HeatExchanger:BadTemps','Temp 1 is greater than Temp 2!'));
    % No heat exchanged
    Q_exchanged = 0;
    dT_LMTD = NaN;
    T_1_out = T_1_in;
    T_2_out = T_2_in;
    T_1 = [];
    T_2 = [];
    Q = [];
    return;
end
if (m_dot_1 <= 0 || m_dot_2 <= 0)
    throw(MException('HeatExchanger:NegativeMassFlow','Negative Massflow in Heat Exchanger'));
end

increments = 20;

if (T_1_in < T_2_in)
    direction = -1;
else
    direction = 1;
end

% Check the phase on incoming fluid
P_crit_1 = CoolProp('Pcrit', '', 0, '', 0, fluid_1);
if (P_1 < P_crit_1)
    T_sat_1 = CoolProp('T', 'P', P_1, 'Q', 1, fluid_1) - 273.15; 
    if (T_sat_1 == T_1_in || T_sat_1 == T_2_in)
        throw(MException('HeatExchanger:TwoPhaseFluid','Fluid 1 enters or leaves two-phase!'));
    end
end
P_crit_2 = CoolProp('Pcrit', '', 0, '', 0, fluid_2);
if (P_2 < P_crit_2)
    T_sat_2 = CoolProp('T', 'P', P_2, 'Q', 1, fluid_2) - 273.15;
    if (T_sat_2 == T_1_in || T_sat_2 == T_2_in)
        throw(MException('HeatExchanger:TwoPhaseFluid','Fluid 2 enters or leaves two-phase!'));
    end
end

h_1_in = CoolProp('HMASS', 'T', T_1_in+273.15, 'P', P_1, fluid_1);
T_1_max = T_2_in;
h_1_max = CoolProp('HMASS', 'T', T_1_max+273.15, 'P', P_1, fluid_1);
T_1_max_practical = T_2_in + direction*dT_pinch;
h_1_max_practical = CoolProp('HMASS', 'T', T_1_max_practical+273.15, 'P', P_1, fluid_1);
h_2_in = CoolProp('HMASS', 'T', T_2_in+273.15, 'P', P_2, fluid_2);
T_2_max = T_1_in;
h_2_max = CoolProp('HMASS', 'T', T_2_max+273.15, 'P', P_2, fluid_2);
T_2_max_practical = T_1_in - direction*dT_pinch;
h_2_max_practical = CoolProp('HMASS', 'T', T_2_max_practical+273.15, 'P', P_2, fluid_2);

Q_1_max = abs( m_dot_1 * (h_1_in - h_1_max) );
Q_2_max = abs( m_dot_2 * (h_2_in - h_2_max) );
Q_1_max_practical = abs( m_dot_1 * (h_1_in - h_1_max_practical) );
Q_2_max_practical = abs( m_dot_2 * (h_2_in - h_2_max_practical) );

if (abs(Q_1_max) < abs(Q_2_max))
    %limitingFluid = 1;
    Q_max = Q_1_max;
    Q_max_practical = Q_1_max_practical;
else
    %limtingFluid = 2;
    Q_max = Q_2_max;
    Q_max_practical = Q_2_max_practical;
end

Q_exchanged = Q_max_practical;
ddT_pinch = 1;
while (ddT_pinch > 0.1)

    dQ = Q_exchanged / increments;
    
    Q(1) = 0;
    h_1(1) = h_1_in;
    T_1(1) = T_1_in;
    h_2(1) = (Q_2_max - Q_exchanged)/m_dot_2 + h_2_max;
    T_2(1) = CoolProp('T', 'HMASS', h_2(1), 'P', P_2, fluid_2) - 273.15;
    dT(1) = direction * (T_1(1) - T_2(1));
    UA(1) = dQ / dT(1);

    for i = 2:increments+1
        Q(i) = Q(i-1) + dQ;
        h_1(i) = h_1(i-1) + dQ/m_dot_1;
        h_2(i) = h_2(i-1) + dQ/m_dot_2;
        T_1(i) = CoolProp('T', 'HMASS', h_1(i), 'P', P_1, fluid_1) - 273.15;
        T_2(i) = CoolProp('T', 'HMASS', h_2(i), 'P', P_2, fluid_2) - 273.15;
        dT(i) = direction * (T_1(i) - T_2(i));
        UA(i) = dQ / dT(i);
    end
    
    min_dT = min(dT);
    ddT_pinch = dT_pinch - min_dT;
    
    % Adjust Q_exchanged
    % Use proportional error approach
    change = ddT_pinch / (T_2_in - T_1_in);
    Q_exchanged = (1-change) * Q_exchanged; 
end

dT_LMTD = Q_exchanged / sum(UA);
effectiveness = Q_exchanged / Q_max;
T_1_out = T_1(end);
T_2_out = T_2(1);

% Check the phase on leaving fluid
if (P_1 < P_crit_1)
    T_sat_1 = CoolProp('T', 'P', P_1, 'Q', 1, fluid_1) - 273.15; 
    if (T_sat_1 > T_1_in && T_sat_1 < T_1_out)
        disp('Caution: Fluid 1 is phase changing in heat exchanger');
    end
    if (T_sat_1 == T_1_out)
        throw(MException('HeatExchanger:TwoPhaseFluid','Fluid 1 leaves two-phase!'));
    end
end
if (P_2 < P_crit_2)
    T_sat_2 = CoolProp('T', 'P', P_2, 'Q', 1, fluid_2) - 273.15;
    if (T_sat_2 > T_2_in && T_sat_2 < T_2_out)
        disp('Caution: Fluid 2 is phase changing in heat exchanger');
    end
    if (T_sat_2 == T_2_out)
        throw(MException('HeatExchanger:TwoPhaseFluid','Fluid 2 leaves two-phase!'));
    end
end

end

