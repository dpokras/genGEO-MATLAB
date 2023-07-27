clear all;
close all;

params = SimulationParameters;
params.coolingMode = 'Wet';
params.orcFluid = 'R245fa';
%params.orcFluid = 'R600a';

% T_in needs higher resolution (0.1) around 150C-200C for each fluid
T_in = 30:1:370;
%T_in = 100:10:370;
%T_in = [66];
T_boil_max_C = MaxSubcritORCBoilTemp(params);
% T_boil_max_C is 123.8 for R245fa
% T_boil_max_C is 107.7 for R600a

Output = zeros(0);
for i = 1:size(T_in, 2)
    T_in_C = T_in(i);
    
    disp(strcat(['Optimizing for input temp: ' num2str(T_in_C, 3)]));
    
    T_boil_C = 40;
    dT_boil = 10;
    dT_boil_old = 10;
    spcc_old = NaN;
    peaks = 0;
    
    dT_pinch = params.dT_orc_pinch;
    dT_dT_pinch = 0;
    dT_dT_pinch_old = -1/0.21;
    peaks_pinch = 0;
    
    while (peaks < 6 && peaks_pinch < 6)
        T_boil_C = T_boil_C + dT_boil;
        dT_pinch = dT_pinch + dT_dT_pinch;
        if (T_boil_C > T_boil_max_C)
            T_boil_C = T_boil_max_C;
        end
        if (T_boil_C < 0)
            % couldn't find any useable values.
            T_boil_C = NaN;
            break;
        end
        
        try
            result_ORC = ORC_Cycle_Tboil(T_in_C, T_boil_C, dT_pinch, params);

            %optimize for 4*100 kg/s
            m_dot = 400;

            Q_preheater_total = m_dot * result_ORC.q_preheater;
            Q_boiler_total = m_dot * result_ORC.q_boiler;
            W_turbine_total = m_dot * result_ORC.w_turbine;
            Q_recuperator_total = 0;
            Q_desuperheater_total = m_dot * result_ORC.q_desuperheater;
            Q_condenser_total = m_dot * result_ORC.q_condenser;
            W_pump_orc_total = m_dot * result_ORC.w_pump;
            W_pump_prod_total = 0;

            W_net_total = W_turbine_total + W_pump_orc_total;

            result_surfacePlant = CapitalCost_SurfacePlant_ORC(Q_preheater_total, Q_boiler_total, W_turbine_total, ...
                Q_recuperator_total, Q_desuperheater_total, Q_condenser_total, ...
                W_pump_orc_total, W_pump_prod_total, ...
                result_ORC.dT_range_CT, result_ORC.dT_LMTD_preheater, result_ORC.dT_LMTD_boiler, 0, params);

            % estimate cost as surface plant + two 40cm diameter wells
            depth = (T_in_C - params.T_surface_air_C)/params.dT_dz;
            C_well = CapitalCost_Well(depth, 0.40, 'Water', 1, params);
            C_total = result_surfacePlant.C_plant + 2*C_well;
            
            result_greenfield = LCOE_Simple(C_total, W_net_total, params);
            spcc = result_greenfield.SpecificCapitalCost;
        catch
            spcc = NaN;
        end
        
        disp(strcat(['Boiling Temp: ' num2str(T_boil_C, '%.1f') ', Pinch: ' num2str(dT_pinch) ', SpCC: ' num2str(spcc*1e3, '%.1f')]));
            
        if (isnan(spcc) && dT_boil > 0)
            % make sure boiling temp is decreasing
            dT_boil = dT_boil * -0.21;
            peaks = peaks + 1;
            disp(strcat(['peak: ' num2str(peaks)]));
        elseif (isnan(spcc) && dT_boil < 0 && ~isnan(spcc_old))
            % it hit an old value that worked, but stopped working again
            % step backwards
            T_boil_C = T_boil_C - dT_boil;
            dT_boil = dT_boil * -0.21;
            peaks = peaks + 1;
            disp(strcat(['peak: ' num2str(peaks)]));
        elseif (isnan(spcc) && dT_dT_pinch ~= 0 && ~isnan(spcc_old))
            dT_pinch = dT_pinch - dT_dT_pinch;
            dT_dT_pinch = dT_dT_pinch * -0.21;
            peaks_pinch = peaks_pinch + 1;
            disp(strcat(['dT_pinch peak: ' num2str(peaks_pinch)]));
        elseif (isnan(spcc) || isnan(spcc_old))
            % do nothing
        elseif (T_boil_C >= T_boil_max_C)
            dT_boil = -0.21 * dT_boil;
            peaks = peaks + 1;
            disp(strcat(['dT_boil peak: ' num2str(peaks)]));
        elseif (dT_pinch < params.dT_orc_pinch)
            dT_dT_pinch = -0.21 * dT_dT_pinch_old;
            dT_dT_pinch_old = dT_dT_pinch;
            dT_pinch = params.dT_orc_pinch - dT_dT_pinch;    
            peaks_pinch = peaks_pinch + 1;
            dT_boil = 0;
            disp(strcat(['dT_pinch peak: ' num2str(peaks_pinch)]));
        elseif (spcc > spcc_old)
            %
            if (peaks >= peaks_pinch)
                % work on pinch variable now
                dT_dT_pinch = -0.21 * dT_dT_pinch_old;
                dT_dT_pinch_old = dT_dT_pinch;
                peaks_pinch = peaks_pinch + 1;
                dT_boil = 0;
                disp(strcat(['dT_pinch peak: ' num2str(peaks_pinch)]));
            else
                % work on boil temp
                dT_boil = -0.21 * dT_boil_old;
                dT_boil_old = dT_boil;
                peaks = peaks + 1;
                dT_dT_pinch = 0;
                disp(strcat(['dT_boil peak: ' num2str(peaks)]));
            end
        end
        
        spcc_old = spcc;
        
    end
    
    disp(strcat(['For input temp: ' num2str(T_in_C, '%.1f') ' found boiling temp of ' num2str(T_boil_C, '%.1f') ' and pinch of ' num2str(dT_pinch)]));
        
    Output = [Output; T_in_C T_boil_C dT_pinch spcc*1e3];
end

subplot(2,1,1);
hold on;
yyaxis left;
plot(Output(:,1),Output(:,2));
ylabel('Boil Temperature [C]');
yyaxis right;
plot(Output(:,1),Output(:,4));
ylabel('SpCC [$/kW]');
hold off;
xlabel('Input Temperature [C]');

subplot(2,1,2);
plot(Output(:,1),Output(:,3));
ylabel('Pinch Temperature [C]');
xlabel('Input Temperature [C]');

% only save first three columns
writematrix(Output(:,1:3), strcat(['data\ORC_Tboil_optimum_minCost_' params.orcFluid '.csv']));