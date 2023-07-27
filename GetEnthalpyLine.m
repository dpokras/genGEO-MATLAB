function [T,s] = GetEnthalpyLine(h_level, T_min, T_max, params)

%h1 = CoolProp('HMASS','P',233e3,'Q',0,'R245fa');
%T_max = CoolProp('T','HMASS',h1,'Q',0,params.orcFluid)-273.15;
dT = (T_max - T_min)/100;

T = [];
s = [];

T_crit = CoolProp('TCRIT', "", 0, "", 0, params.orcFluid)-273.15;

for T_val = T_min:dT:T_max
    % Gotta do it this funny way because coolprop won't take in T,H to give
    % s. So we have to look up many s and then interpolate.
    if (T_val < T_crit)
        s_min = CoolProp('SMASS','T',T_val+273.15,'Q',0,params.orcFluid);
    else
        s_min = 1.6e3;
    end
    s_max = 1.9e3;
    ds = (s_max - s_min) / 100;
    s_range = s_min:ds:s_max;
    for i = 1:size(s_range,2)
        h_range(i) = CoolProp('HMASS', 'T', T_val+273.15, 'SMASS', s_range(i), params.orcFluid);
    end
    s_val = interp1(h_range,s_range,h_level);
    T = [T T_val];
    s = [s s_val];
end


end