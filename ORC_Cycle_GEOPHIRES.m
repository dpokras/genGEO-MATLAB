function result = ORC_Cycle_GEOPHIRES(T_in_C, params)
   
    type = 'SubcriticalORC';
    %type = 'SupercriticalORC';
    %type = 'SingleFlash';
    %type = 'DoubleFlash';
    
    T_ambient_C = params.T_surface_air_C;
    
    if (strcmp(type,'SubcriticalORC'))
        if (T_ambient_C < 15)
            C1 = 2.746E-3;
            C0 = -8.3806E-2;
            D1 = 2.713E-3;
            D0 = -9.1841E-2;
            Tfraction = (T_ambient_C-5)/10;
        else
            C1 = 2.713E-3;
            C0 = -9.1841E-2;
            D1 = 2.676E-3;
            D0 = -1.012E-1;
            Tfraction = (T_ambient_C-15)/10;
        end
        etaull = C1*TenteringPP + C0;
        etauul = D1*TenteringPP + D0;
        etau = (1-Tfraction)*etaull + Tfraction*etauul;
        if (T_ambient_C < 15)
            C1 = 0.0894;
            C0 = 55.6;
            D1 = 0.0894;
            D0 = 62.6;
            Tfraction = (T_ambient_C-5)/10;
        else
            C1 = 0.0894;
            C0 = 62.6;
            D1 = 0.0894;
            D0 = 69.6;
            Tfraction = (T_ambient_C-15)/10;
        end
        reinjtll = C1*TenteringPP + C0;
        reinjtul = D1*TenteringPP + D0;
        ReinjTemp = (1-Tfraction)*reinjtll + Tfraction*reinjtul;
    elseif (strcmp(type,'SupercriticalORC'))
        if (T_ambient_C < 15)
            C2 = -1.55E-5;    
            C1 = 7.604E-3;
            C0 = -3.78E-1;
            D2 = -1.499E-5;
            D1 = 7.4268E-3;
            D0 = -3.7915E-1;
            Tfraction = (T_ambient_C-5)/10;
        else
            C2 = -1.499E-5;        
            C1 = 7.4268E-3;
            C0 = -3.7915E-1;
            D2 = -1.55E-5;
            D1 = 7.55136E-3;
            D0 = -4.041E-1;
            Tfraction = (T_ambient_C-15)/10;
        end
        etaull = C2*TenteringPP^2 + C1*TenteringPP + C0;
        etauul = D2*TenteringPP^2 + D1*TenteringPP + D0;
        etau = (1-Tfraction)*etaull + Tfraction*etauul;
        if (T_ambient_C < 15)
            C1 = 0.02;
            C0 = 49.26;
            D1 = 0.02;
            D0 = 56.26;
            Tfraction = (T_ambient_C-5)/10;
        else
            C1 = 0.02;
            C0 = 56.26;
            D1 = 0.02;
            D0 = 63.26;
            Tfraction = (T_ambient_C-15)/10;
        end
        reinjtll = C1*TenteringPP + C0;
        reinjtul = D1*TenteringPP + D0;
        ReinjTemp = (1-Tfraction)*reinjtll + Tfraction*reinjtul;
    elseif (strcmp(type,'SingleFlash'))
        if (T_ambient_C < 15)
            C2 = -4.27318E-7;
            C1 = 8.65629E-4;
            C0 = 1.78931E-1;
            D2 = -5.85412E-7;
            D1 = 9.68352E-4;
            D0 = 1.58056E-1;
            Tfraction = (T_ambient_C-5)/10;
        else
            C2 = -5.85412E-7;
            C1 = 9.68352E-4;
            C0 = 1.58056E-1;
            D2 = -7.78996E-7;
            D1 = 1.09230E-3;
            D0 = 1.33708E-1;
            Tfraction = (T_ambient_C-15)/10;
        end
        etaull = C2*TenteringPP^2 + C1*TenteringPP + C0;
        etauul = D2*TenteringPP^2 + D1*TenteringPP + D0;
        etau = (1-Tfraction)*etaull + Tfraction*etauul;
        if (T_ambient_C < 15)
            C2 = -1.11519E-3;
            C1 = 7.79126E-1;
            C0 = -10.2242;
            D2 = -1.10232E-3;
            D1 = 7.83893E-1;
            D0 = -5.17039;
            Tfraction = (T_ambient_C-5)/10;
        else
            C2 = -1.10232E-3;
            C1 = 7.83893E-1;
            C0 = -5.17039;
            D2 = -1.08914E-3;
            D1 = 7.88562E-1;
            D0 = -1.89707E-1;
            Tfraction = (T_ambient_C-15)/10;
        end
        reinjtll = C2*TenteringPP^2 + C1*TenteringPP + C0;
        reinjtul = D2*TenteringPP^2 + D1*TenteringPP + D0;
        ReinjTemp = (1-Tfraction)*reinjtll + Tfraction*reinjtul;
    elseif (strcmp(type,'DoubleFlash'))
        if (T_ambient_C < 15)
            C2 = -1.200E-6;
            C1 = 1.22731E-3;
            C0 = 2.26956E-1;
            D2 = -1.42165E-6;
            D1 = 1.37050E-3;
            D0 = 1.99847E-1;
            Tfraction = (T_ambient_C-5)/10;
        else
            C2 = -1.42165E-6;
            C1 = 1.37050E-3;
            C0 = 1.99847E-1;
            D2 = -1.66771E-6;
            D1 = 1.53079E-3;
            D0 = 1.69439E-1;
            Tfraction = (T_ambient_C-15)/10;
        end
        etaull = C2*TenteringPP^2 + C1*TenteringPP + C0;
        etauul = D2*TenteringPP^2 + D1*TenteringPP + D0;
        etau = (1-Tfraction)*etaull + Tfraction*etauul;
        if (T_ambient_C < 15)
            C2 = -7.70928E-4;
            C1 = 5.02466E-1;
            C0 = 5.22091;
            D2 = -7.69455E-4;
            D1 = 5.09406E-1;
            D0 = 11.6859;
            Tfraction = (T_ambient_C-5)/10;
        else
            C2 = -7.69455E-4;
            C1 = 5.09406E-1;
            C0 = 11.6859;
            D2 = -7.67751E-4;
            D1 = 5.16356E-1;
            D0 = 18.0798;
            Tfraction = (T_ambient_C-15)/10;
        end
        reinjtll = C2*TenteringPP^2 + C1*TenteringPP + C0;
        reinjtul = D2*TenteringPP^2 + D1*TenteringPP + D0;
        ReinjTemp = (1-Tfraction)*reinjtll + Tfraction*reinjtul;
    else
        throw(MException('ORC_Cycle_GEOPHIRES:UnknownOrcModel','Unknown GEOPHIRES ORC Model'));
    end

    % Exergy At surface (instead of using geophires correlation)
    h_0 = CoolProp('HMASS','P',P_inj_surface,'T',T_ambient_C+273.15,fluid);
    s_0 = CoolProp('SMASS','P',P_inj_surface,'T',T_ambient_C+273.15,fluid);
    s_prod_surface = CoolProp('SMASS','P',P_prod_surface,'T',T_in_C+273.15,fluid);
    % Available exergy relationship calculation here is nearly identical to
    % regression relationship in geophires.
    Available_Exergy = h_prod_surface - h_0 - (T_ambient_C+273.15)*(s_prod_surface - s_0);

    P_sat = CoolProp('P', 'T', T_in_C + 273.15, 'Q', 0, 'Water');
    cp = CoolProp('CPMASS', 'T', T_in_C + 273.15, 'P', P_sat + 100e3, 'Water');
    q_boiler = cp * (T_in_C - T_ambient_C);
    w_turbine = Available_Exergy * etau;
    q_condenser = -1* (q_boiler - w_turbine);
    w_net = w_turbine;
    
    T_out_C = ReinjTemp;
    
    result.q_boiler = q_boiler;
    result.w_turbine = w_turbine;
    result.q_condenser = q_condenser;
    result.w_net = w_net;
    result.T_out_C = T_out_C;
    
end