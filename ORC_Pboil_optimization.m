T_ambient_C = 15;
dT_approach = 7;
dT_pinch = 5;
eta_pump = 0.9;
eta_turbine = 0.8;
coolingMode = 'Wet';

T_in = 170:1:370;
P_crit = CoolProp('PCRIT', "", 0, "", 0, 'R245fa');

Output = zeros(0);
for i = 1:size(T_in, 2)
    T_in_C = T_in(i);
    
    disp(strcat(['Optimizing for input temp: ' num2str(T_in_C, 3)]));
    
    P_boil_Pa = 4e6;
    dP_boil_Pa = 1e5;
    dw_net = 100;
    w_net_old = 0;
    peaks = 0;
    
    while (peaks < 6)
        P_boil_Pa = P_boil_Pa + dP_boil_Pa;
        if (P_boil_Pa < P_crit + 1e3)
            P_boil_Pa = P_crit + 1e3;
            disp('setting boil to Pcrit');
        end
        
        try
            disp(strcat(['Trying a boiling pressure of ' num2str(P_boil_Pa/1e6,'%.3f') ' MPa.']));
            [q_boiler, w_turbine, q_recuperator, q_desuperheater, q_condenser, w_pump, w_cooler, w_condenser, w_net, ...
                dT_range_CT, dT_LMTD_boiler, dT_LMTD_recuperator, T_out_C] ...
                = ORC_Cycle_Supercrit_Pboil(T_in_C, P_boil_Pa, T_ambient_C, dT_approach, dT_pinch, eta_pump, eta_turbine, coolingMode);
            disp(strcat(['Boiling Pressure: ' num2str(P_boil_Pa/1e6, '%.3f') ' MPa, w_net: ' num2str(w_net/1e3, '%.3f') ' kWe.']));
            
            if ((isnan(w_net_old)==false && w_net < w_net_old) || P_boil_Pa == P_crit + 1e3)
                dP_boil_Pa = dP_boil_Pa * -0.21;
                peaks = peaks + 1;
                disp(strcat(['peak: ' num2str(peaks)]));
            end
        
            w_net_old = w_net;
        catch ME
            if (strcmp(ME.identifier,'ORC_Cycle_Supercrit_Pboil:lowInletTemp'))
                % Make sure P_boil is decreasing
                if (dP_boil_Pa > 0)
                    dP_boil_Pa = dP_boil_Pa * -0.21;
                    peaks = peaks + 1;
                    disp(strcat(['peak: ' num2str(peaks)]));
                end
                
                w_net_old = NaN;
            end
        end
            

        
    end
    
    disp(strcat(['For input temp: ' num2str(T_in_C, '%.3f') ' found boiling pressure of ' num2str(P_boil_Pa/1e3, '%.3f')]));
        
    Output = [Output; T_in_C P_boil_Pa];
    
    writematrix(Output, 'data\ORC_Pboil_optimum.csv');
end