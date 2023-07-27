function [f_cooling, f_condensing] = ParasiticPowerFraction_CoolingTower(T_ambient_C, dT_approach_CT, dT_range_CT, coolingMode)

    if (strcmp(coolingMode,'Wet'))
        %wet
        a_cool_wet = 1.20e0;
        b_cool_wet = 0;
        c_cool_wet = -3.79e-3;
        d_cool_wet = 1.95e-2;
        f_cooling = a_cool_wet*(1/dT_approach_CT) + b_cool_wet*(T_ambient_C+273.15) + c_cool_wet*(T_ambient_C+273.15)/dT_approach_CT + d_cool_wet*(1/(dT_approach_CT + dT_range_CT)); 
        
        a_cond_wet = 1.65e0;
        b_cond_wet = -6.24e-6;
        c_cond_wet = -5.03e-3;
        d_cond_wet = 0;
        f_condensing = a_cond_wet*(1/dT_approach_CT) + b_cond_wet*(T_ambient_C+273.15) + c_cond_wet*(T_ambient_C+273.15)/dT_approach_CT + d_cond_wet*(1/(dT_approach_CT + dT_range_CT)); 
    
    elseif (strcmp(coolingMode,'Dry'))
        %dry
        a_cool_dry = 7.65e-1;
        b_cool_dry = 0;
        c_cool_dry = 0;
        d_cool_dry = 1.28e-1;
        f_cooling = a_cool_dry*(1/dT_approach_CT) + b_cool_dry*(T_ambient_C+273.15) + c_cool_dry*(T_ambient_C+273.15)/dT_approach_CT + d_cool_dry*(1/(dT_approach_CT + dT_range_CT)); 
       
        a_cond_dry = 6.19e-1;
        b_cond_dry = 0;
        c_cond_dry = 0;
        d_cond_dry = 0;
        f_condensing = a_cond_dry*(1/dT_approach_CT) + b_cond_dry*(T_ambient_C+273.15) + c_cond_dry*(T_ambient_C+273.15)/dT_approach_CT + d_cond_dry*(1/(dT_approach_CT + dT_range_CT)); 
    else
        throw(MException('ParasiticPowerFraction:UnknownCoolingMode','Unknown Cooling Mode'));
    end
end

