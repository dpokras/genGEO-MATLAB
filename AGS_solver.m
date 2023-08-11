function [result,x_opt] = AGS_solver(useMySolver, x0, params)
    
    time = [0.1,0.5,1:1:30];
    time_total = 3600 * 24 * 365 * time(end); %s
    Energy_net = [];
    Thermal_Energy_net = [];
    Energy_reservoir = [];
    Energy_sold = [];
    Power = [];
    
    params.time_years = time(end); % set to first time value or any other suitable value
    
    if useMySolver % compute x_opt once outside of the loop
        if params.config == 1 || params.config == 2
            myCostFunc = @(x)costfunc_config1(total_analytic_system_co2(wrapMyParameters(params,x)));
        elseif params.config || params.config == 3 || params.config == 4
            myCostFunc = @(x)costfunc_config3(total_analytic_system_co2(wrapMyParameters(params,x)));
        end

        myNonLinConstraints = @(x)contrfunc(total_analytic_system_co2(wrapMyParameters(params,x)));

        myOptions = optimoptions('fmincon',...
            'Algorithm','interior-point',...
            'ConstraintTolerance',1.0000e-03,...
            'StepTolerance',1.0000e-03,...
            'MaxIterations',30);

        x_opt = fmincon( ...
            myCostFunc, ...
            x0, ...
            [],[], ...
            [], [], ...
           [10],[inf],myNonLinConstraints,myOptions);
    else
        x_opt = [x0]; % for useMySolver = 0, use initial guess
    end

    for i = 1:size(time,2)
        params.time_years = time(i);
        params.m_dot = x_opt(1); % use optimized m_dot

        if i == 1
            time_seconds = 3600 * 24 * 365 * time(i); %s
        else
            time_seconds = 3600 * 24 * 365 * (time(i) - time(i-1)); %s
        end

        result = total_analytic_system_co2(wrapMyParameters(params,x_opt));

        Energy_net(end + 1) = result.W_net .* time_seconds;  %J
        Thermal_Energy_net(end + 1) = result.Q_net .* time_seconds;  %J
        Energy_reservoir(end + 1) = result.delta_H(end) .* time_seconds;  %J
        Energy_sold(end + 1) = (result.W_net + result.Q_net) .* time_seconds;  %J
        Power(end + 1) = (result.W_net + result.Q_net);  %W
    end

    result.opt_m_dot = x_opt(1);
    result.S_ratio = result.S_ratio;
    result.Power = Power;  %W
    result.Power_electric_sold_avg = sum(Energy_net)/time_total;
    result.Power_reservoir_avg = sum(Energy_reservoir)/time_total;
    result.Power_heat_sold_avg = sum(Thermal_Energy_net)/time_total;
    result.Power_total_sold_avg = sum(Energy_sold)/time_total;
    disp(['Net Electric Power (Time Averaged): [' num2str(result.Power_electric_sold_avg/1e3) '] kW']);
    disp(['Net Thermal Power (Time Averaged): [' num2str(result.Power_heat_sold_avg/1e3) '] kW']);
    disp(['Net Power Sold (Time Averaged): [' num2str(result.Power_total_sold_avg/1e3) '] kW']);   
    disp(['Net Electic Power: [' num2str(result.W_net/1e3) '] kW']);
    disp(['Net Thermal Power: [' num2str(result.Q_net/1e3) '] kW']);
    disp(['Net Power Reservoir (Time Averaged): [' num2str(result.Power_reservoir_avg/1e3) '] kW']);
    disp(['SpCC Power Brownfield: [' num2str(result.SpCC_W_brownfield) '] 2022$/W']);
    disp(['SpCC Heat Brownfield: [' num2str(result.SpCC_Q_brownfield) '] 2022$/W']);
    disp(['SpCC dH Brownfield: [' num2str(result.SpCC_dH_greenfield) '] 2022$/W']);
    disp(['SpCC Power Greenfield: [' num2str(result.SpCC_W_greenfield) '] 2022$/W']);
    disp(['SpCC Heat Greenfield: [' num2str(result.SpCC_Q_greenfield) '] 2022$/W']);
    disp(['SpCC dH Greenfield: [' num2str(result.SpCC_dH_greenfield) '] 2022$/W']);
    disp(['Depth: [' num2str(params.depth) '] m']);
    disp(['Lateral Length: [' num2str(params.res_length) '] m']);
    disp(['no of side streams: [' num2str(params.n_streams) ']']);
    disp(['Mass Flow: [' num2str(x_opt(1)) '] kg/s']);
    disp(['max speed: [' num2str(result.max_speed) '] m/s']);
    disp(['entopy into turbine: [' num2str(result.s_turb_in) '] J/kgK']);
end