T_ambient_C = 15;
dT_approach = 7;
dT_pinch = 5;
eta_pump = 0.9;
eta_turbine = 0.8;
coolingMode = 'Wet';

T_in = 30:1:370;

%Output = zeros(0);
Output = [{'T_in_C'} {'T_ambient_C'} {'dT_approach'} {'dT_pinch'} {'eta_pump'} {'eta_turbine'} {'coolingMode'} ...
        {'q_preheater'} {'q_boiler'} {'w_turbine'} {'q_recuperator'} {'q_desuperheater'} {'q_condenser'} {'w_pump'} {'w_cooler'} {'w_condenser'} {'w_net'} ...
        {'dT_range_CT'} {'dT_LMTD_preheater'} {'dT_LMTD_boiler'} {'dT_LMTD_recuperator'} {'T_out_C'}];

    
for i = 1:size(T_in, 2)
    T_in_C = T_in(i);
    
    disp(strcat(['Running for input temp: ' num2str(T_in_C, 3)]));
    
    [q_preheater, q_boiler, w_turbine, q_recuperator, q_desuperheater, q_condenser, w_pump, w_cooler, w_condenser, w_net, ...
        dT_range_CT, dT_LMTD_preheater, dT_LMTD_boiler, dT_LMTD_recuperator, T_out_C] ...
        = ORC_Cycle(T_in_C, T_ambient_C, dT_approach, dT_pinch, eta_pump, eta_turbine, coolingMode);
    
    disp(strcat(['For input temp: ' num2str(T_in_C, '%.1f') ' found T_out of ' num2str(T_out_C, '%.1f')]));
     
    Output = [Output; T_in_C T_ambient_C dT_approach dT_pinch eta_pump eta_turbine {coolingMode} ...
        q_preheater q_boiler w_turbine q_recuperator q_desuperheater q_condenser w_pump w_cooler w_condenser w_net ...
        dT_range_CT dT_LMTD_preheater dT_LMTD_boiler dT_LMTD_recuperator T_out_C];
end

xlswrite('ORC_Lookup.csv', Output);

%clear persistent variable
%clear ORC_Cycle;