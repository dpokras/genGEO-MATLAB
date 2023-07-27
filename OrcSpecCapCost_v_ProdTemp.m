%need to clear all to remove persistent variables in ORC_Cycle
clear all;
close all;

params = SimulationParameters;
%params.orcFluid = 'R245fa';
params.orcFluid = 'R600a';


prodTemp = 80:1:320;%320;
capacity = [5, 10, 25, 50, 75, 100];

model = 'PROPOSED';
%model = 'GEOPHIRES';


%GEOPHIRES Cycle
cycleType = 'SubcriticalORC';
%cycleType = 'SupercriticalORC';
%cycleType = 'SingleFlash';
%cycleType = 'DoubleFlash';

if (strcmp(model,'PROPOSED'))
    cycleType = 'SubcriticalORC';
end

Output = zeros(0);
for i = 1:size(prodTemp,2)
   
    T_in_C = prodTemp(i);
    result = ORC_Cycle(T_in_C, params);
    
    for j = 1:size(capacity, 2)
        
        disp(strcat(['Iteration x=' num2str(prodTemp(i)) ', y=' num2str(capacity(j)) '.']));
        
        cap_MW = capacity(j);
        % find m_dot to meet capacity
        m_dot = cap_MW * 1e6 / result.w_net;
        
        Q_preheater = m_dot * result.q_preheater;
        Q_boiler = m_dot * result.q_boiler;
        W_turbine = m_dot * result.w_turbine;
        Q_recuperator = m_dot * result.q_recuperator;
        Q_desuperheater = m_dot * result.q_desuperheater;
        Q_condenser = m_dot * result.q_condenser;
        W_pump_orc = m_dot * result.w_pump;
        W_cooler_orc = m_dot * result.w_cooler;
        W_condenser_orc = m_dot * result.w_condenser;
        W_pump_prod = 0;
        
        if (strcmp(model,'PROPOSED'))    
            result_cost = CapitalCost_SurfacePlant_ORC(Q_preheater, Q_boiler, W_turbine, ...
                Q_recuperator, Q_desuperheater, Q_condenser, ...
                W_pump_orc, W_pump_prod, result.dT_range_CT, result.dT_LMTD_preheater, ...
                result.dT_LMTD_boiler, 0, params);
            C_surfacePlant = result_cost.C_plant;
        elseif (strcmp(model,'GEOPHIRES'))
            ElectricityProduced_MW = cap_MW;
            C_plant_geophires = CapitalCost_SurfacePlant_ORC_GEOPHIRES(T_in_C, ElectricityProduced_MW, cycleType, params);
            C_surfacePlant = C_plant_geophires;
        end
        
        SpecificCapitalCost = C_surfacePlant / (cap_MW * 1000);

        Output = [ Output; T_in_C cap_MW {model} {cycleType} SpecificCapitalCost ];
    end
end

hold on;
for i = 1:size(capacity,2)
    cap_MW = capacity(i);
    ind = cell2mat(Output(:,2)) == cap_MW;
    plot(cell2mat(Output(ind,1)),cell2mat(Output(ind,5)),'DisplayName',strcat([ num2str(cap_MW) ' MWe']));
end
grid on;
ylabel('Specific Power [$/kWe]');
xlabel('Input Temperature [\circC]');
legend();

filename = strcat(['data\SpecCapCost_' datestr(now,'yyyymmdd_HHMMSS') '.xlsx']);
xlswrite(filename, Output);
