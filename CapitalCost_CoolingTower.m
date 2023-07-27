function C_coolingTowers = CapitalCost_CoolingTower(Q_cooler, Q_condenser, TDC, dT_range_CT, params)

% Heat in Watts

if (strcmp(params.coolingMode,'Wet'))
    a_cool = 5.58e3;
    b_cool = 0;
    c_cool = -1.77e1;
    d_cool = 1.96e2;
    c_cooling = a_cool*(1/params.dT_approach) + b_cool*(params.T_surface_air_C+273.15) + c_cool*(params.T_surface_air_C+273.15)/params.dT_approach + d_cool*(1/(params.dT_approach + dT_range_CT)); 
    a_cond = 4.08e3;
    b_cond = -1.54e-2;
    c_cond = -1.24e1;
    d_cond = 0;
    c_condensing = a_cond*(1/params.dT_approach) + b_cond*(params.T_surface_air_C+273.15) + c_cond*(params.T_surface_air_C+273.15)/params.dT_approach + d_cond*(1/(params.dT_approach + dT_range_CT)); 
elseif (strcmp(params.coolingMode,'Dry'))
    a_cool = 7.31e3;
    b_cool = 0;
    c_cool = 0;
    d_cool = 1.23e3;
    c_cooling = a_cool*(1/params.dT_approach) + b_cool*(params.T_surface_air_C+273.15) + c_cool*(params.T_surface_air_C+273.15)/params.dT_approach + d_cool*(1/(params.dT_approach + dT_range_CT)); 
    a_cond = 1.91e3;
    b_cond = 0;
    c_cond = 0;
    d_cond = 0;
    c_condensing = a_cond*(1/params.dT_approach) + b_cond*(params.T_surface_air_C+273.15) + c_cond*(params.T_surface_air_C+273.15)/params.dT_approach + d_cond*(1/(params.dT_approach + dT_range_CT)); 
else
    throw(MException('CapitalCost_CoolingTower:UnknownCoolingMode','Unknown Cooling Mode'));
end
    
    
% c_cooling_wet and c_condensing_wet both in units of $/kWth

% Reference case 1000 kWth (1e6 Wth)
Q_Ref_BAC = 1e6;
F_cooling = abs(Q_cooler)/(abs(Q_cooler)+abs(Q_condenser));
C_Ref_BAC = Q_Ref_BAC * TDC * (F_cooling*(c_cooling/1e3) + (1-F_cooling)*(c_condensing/1e3));
C_coolingTowers = C_Ref_BAC * (abs(Q_cooler+Q_condenser)/Q_Ref_BAC)^0.8;

end

