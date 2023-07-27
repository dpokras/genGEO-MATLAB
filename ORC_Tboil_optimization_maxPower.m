clear all;
close all;

params = SimulationParameters;
params.coolingMode = 'Wet';
%params.orcFluid = 'R245fa';
params.orcFluid = 'R600a';

T_in = 30:1:370;
T_boil_max_C = MaxSubcritORCBoilTemp(params);
% T_boil_max_C is 123.8 for R245fa
% T_boil_max_C is 107.7 for R600a

dT_pinch = params.dT_orc_pinch;

Output = zeros(0);
for i = 1:size(T_in, 2)
    T_in_C = T_in(i);
    
    disp(strcat(['Optimizing for input temp: ' num2str(T_in_C, 3)]));
    
    T_boil_C = 40;
    dT_boil = 10;
    dw_net = 100;
    w_net_old = 0;
    peaks = 0;
    
    while (peaks < 6)
        T_boil_C = T_boil_C + dT_boil;
        if (T_boil_C > T_boil_max_C)
            T_boil_C = T_boil_max_C;
        end
        
        try
            result = ORC_Cycle_Tboil(T_in_C, T_boil_C, dT_pinch, params);
            w_net = result.w_net;
        catch
            w_net = NaN;
        end
        
        disp(strcat(['Boiling Temp: ' num2str(T_boil_C, '%.1f') ' , w_net: ' num2str(w_net, '%.1f')]));
            
        if (isnan(w_net) && dT_boil > 0)
            % make sure boiling temp is decreasing
            dT_boil = dT_boil * -0.21;
            peaks = peaks + 1;
            disp(strcat(['peak: ' num2str(peaks)]));
        elseif (w_net < w_net_old || T_boil_C == T_boil_max_C)
            dT_boil = dT_boil * -0.21;
            peaks = peaks + 1;
            disp(strcat(['peak: ' num2str(peaks)]));
        end
        
        w_net_old = w_net;
        
    end
    
    disp(strcat(['For input temp: ' num2str(T_in_C, '%.1f') ' found boiling temp of ' num2str(T_boil_C, '%.1f')]));
        
    Output = [Output; T_in_C T_boil_C dT_pinch];
end

writematrix(Output, strcat(['data\ORC_Tboil_optimum_maxPower_' params.orcFluid '.csv']));