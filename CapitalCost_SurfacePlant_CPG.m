function C_SurfacePlant = CapitalCost_SurfacePlant_CPG(result_turbine, result_HEX, result_pumps, result_cooling_tower, result_tank, result_comp_unit, params)

%% Calculate all the heat exchanger costs

C_BM_hex = CapitalCost_SurfacePlantfunc.HEX([result_HEX.A_1, result_HEX.A_2], [result_HEX.P_tubeside_1, result_HEX.P_tubeside_2],params);
C_BM_hex_comp = CapitalCost_SurfacePlantfunc.HEX(result_comp_unit.A, result_comp_unit.P_tubeside, params);

%% Calcualte the turbine costs
if result_turbine.W_turbine > 0
    C_BM_turbine = CapitalCost_SurfacePlantfunc.turbine(result_turbine, params);
else
    C_BM_turbine = 0;
end

%% Calculate the pumping costs
if result_pumps.result_cooling_pump.m_dot > 0
    C_BM_cooling_water_pump = CapitalCost_SurfacePlantfunc.pump(result_pumps.result_cooling_pump, params);
else
    C_BM_cooling_water_pump = 0;
end

%% Calculate the cooling tower costs
if result_cooling_tower.W_cooling_tower > 0
    C_BM_cooling_tower = CapitalCost_SurfacePlantfunc.cooling_tower(result_cooling_tower, params);
else
    C_BM_cooling_tower = 0;
end

%% Auxillary and start_up costs
%% Calculate the compressor costs

C_BM_compressor = CapitalCost_SurfacePlantfunc.compressor(result_comp_unit, params);

%% Calculate pump costs
%% Calculate the pumping costs
if result_pumps.result_SU_pump.m_dot > 0
    C_BM_SU_pump = CapitalCost_SurfacePlantfunc.pump(result_pumps.result_SU_pump, params);
else
    C_BM_SU_pump = 0;
end

%% Calculate the pumping costs
if result_pumps.result_filling_pump.m_dot > 0
    C_BM_filling_pump = CapitalCost_SurfacePlantfunc.pump(result_pumps.result_filling_pump, params);
else
    C_BM_filling_pump = 0;
end

%% Calculate the storage tanks

C_BM_storage_tank = CapitalCost_SurfacePlantfunc.tank(result_tank, params);

%% Calulate total bare-module costs

C_TBM = sum(C_BM_hex) + C_BM_turbine + C_BM_cooling_water_pump + C_BM_cooling_tower + C_BM_storage_tank + C_BM_SU_pump + C_BM_filling_pump; %USD

C_SurfacePlant = CapitalCost_SurfacePlantfunc.external_costs(C_TBM); %USD
end

