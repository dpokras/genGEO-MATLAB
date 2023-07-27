function [q_exchanged_orc, q_exchanged_geo, dT_LMTD, T_orc_out, T_geo_out, T_orc, T_geo, Q, m_dot_orc, m_dot_geo, mdot_ratio] ...
    = HeatExchanger_OptMdot(T_orc_in, P_orc, fluid_orc, T_geo_in, P_geo, fluid_geo, dT_pinch, T_orc_out_min)


%maximizeHeatFromStream = 'orc';
maximizeHeatFromStream = 'geo';


mdot_ratio = 1;
d_mdot_ratio = 0.2;
q_exchanged_old = NaN;
peaks = 0;

while (peaks < 4)

    mdot_ratio = mdot_ratio + d_mdot_ratio;
    
    % mdot_ratio = m_dot_geo / m_dot_orc;
    m_dot_orc = 1;
    m_dot_geo = m_dot_orc * mdot_ratio;
    %disp(strcat(['Trying mdot_ratio ' num2str(mdot_ratio) ', m_dot_2 ' num2str(m_dot_geo)]));
        
    [Q_exchanged, dT_LMTD, T_orc_out, T_geo_out, T_orc, T_geo, Q] = HeatExchanger(T_orc_in, P_orc, m_dot_orc, fluid_orc, T_geo_in, P_geo, m_dot_geo, fluid_geo, dT_pinch);
    
    q_exchanged_orc = Q_exchanged / m_dot_orc;
    q_exchanged_geo = Q_exchanged / m_dot_geo;
    
    %disp(strcat(['Obtained q_exchanged_1 ' num2str(q_exchanged_orc) ', q_exchanged_2 ' num2str(q_exchanged_geo)]));
    
    if (strcmp(maximizeHeatFromStream, 'orc'))
        q_exchanged = q_exchanged_orc;
    elseif (strcmp(maximizeHeatFromStream, 'geo'))
        q_exchanged = q_exchanged_geo;
    else
        throw(MException('HeatExchanger_OptMdot:unknownOptimizeMethod','Unknown Optimize Method'));
    end
    
    %if temp too low, make sure mdot_ratio is increasing
    if (~isnan(T_orc_out_min) && T_orc_out < T_orc_out_min && d_mdot_ratio < 0)
        d_mdot_ratio = -0.21 * d_mdot_ratio;
        peaks = peaks + 1;
        %disp('peaks + 1');
        q_exchanged = NaN;
    elseif (~isnan(q_exchanged) && ~isnan(q_exchanged_old) && q_exchanged < q_exchanged_old)
        d_mdot_ratio = -0.21 * d_mdot_ratio;
        peaks = peaks + 1;
        %disp('peaks + 1');
    end
    
    q_exchanged_old = q_exchanged;
end


end

