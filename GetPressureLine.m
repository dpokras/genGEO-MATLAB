function [T,s] = GetPressureLine(P, T_min, T_max, params)

T_sat = CoolProp('T','P',P,'Q',0,params.orcFluid)-273.15;
dT_lower = (T_sat - T_min)/100;
dT_upper = (T_max - T_sat)/100;

T = [];
s = [];
% lower
% Don't plot lower lines, they are just on the dome
% for T_lower=T_min:dT_lower:(T_sat-dT_lower)
%     s_lower = CoolProp('SMASS', 'T', T_lower+273.15, 'P', P, params.orcFluid);
%     T = [T T_lower];
%     s = [s s_lower];
% end
s_lower = CoolProp('SMASS', 'T', T_sat+273.15, 'Q', 0, params.orcFluid);
T = [T T_sat];
s = [s s_lower];
% upper
s_upper = CoolProp('SMASS', 'T', T_sat+273.15, 'Q', 1, params.orcFluid);
T = [T T_sat];
s = [s s_upper];
for T_upper=(T_sat+dT_upper):dT_upper:T_max
    s_upper = CoolProp('SMASS', 'T', T_upper+273.15, 'P', P, params.orcFluid);
    T = [T T_upper];
    s = [s s_upper];
end
end