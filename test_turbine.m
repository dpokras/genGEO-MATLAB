params = SimulationParameters;

P_turb_in = 17e6
T_turb_in = 120


%%% States
%%State 1
h_turb_in = CoolProp.PropsSI('HMASS', 'P', P_turb_in, 'T', T_turb_in+273.15, params.fluid);
s_turb_in = CoolProp.PropsSI('SMASS', 'P', P_turb_in, 'T', T_turb_in+273.15, params.fluid);
disp(h_turb_in/1e3)

%%State 2
P_turb_out_iter = 74e5;
toggle = 0;
go = 1;
iter = 1;

while go == 1

%                 disp(strcat(['interation: ' num2str(iter)]));
    
    %isentropic output
    h_turb_out_isen = CoolProp.PropsSI('HMASS', 'P', P_turb_out_iter, 'SMASS', s_turb_in, params.fluid);
    %efficiency losses
    h_turb_out_irr = h_turb_in - params.eta_cpg_turbine * (h_turb_in - h_turb_out_isen);
    s_turb_out_irr = CoolProp.PropsSI('SMASS', 'P', P_turb_out_iter, 'HMASS', h_turb_out_irr, params.fluid);
    T_turb_out_irr = CoolProp.PropsSI('T', 'P', P_turb_out_iter, 'SMASS', s_turb_out_irr, params.fluid);
    phase = CoolProp.PhaseSI('P', P_turb_out_iter, 'SMASS', s_turb_in, params.fluid)
    
%                 disp(strcat(['Turbine Temperature Out: ' num2str(T_turb_out_irr-273.15) 'C']));
%                 disp(strcat(['Turbine Pressure Out: ' num2str(P_turb_out_iter/1e5) ' bar']));
%                 disp(phase);

    P_turb_out_iter = P_turb_out_iter - 0.1e5
    iter = iter + 1

    % Turbine outlet pressure is determined by the while loop above,
    % and a 0.5 bar buffer is added to avoid phase change and 
    % damage to the machinery.
    %
    % PhaseChangeAllowed toggle means if the turbine can handle
    % phase change from supercritical fluid to liquid.

    if toggle == 1
        go = 0;
    elseif params.PhaseChangeAllowed == false
        if strcmp(string(phase), 'liquid') == 1 || contains(string(phase), 'twophase') == 1

            P_turb_out_iter = P_turb_out_iter + 0.6e5;
            toggle = 1;
        end
    elseif params.PhaseChangeAllowed == true
        if contains(string(phase), 'twophase') == 1
            
            P_turb_out_iter = P_turb_out_iter + 0.6e5;
            toggle = 1;
        end                       
    end
end
P_turb_out_irr = P_turb_out_iter;

try
W_turbine = params.m_dot * (h_turb_in - h_turb_out_irr);
catch ME
    if (W_turbine < 0)
        causeException = (MException('total_analytic_system_CO2:TurbinePowerNegative','Turbine Power is Negative'));
        ME = addCause(ME,causeException);
    end
end