function result = CapitalCost_SurfacePlant_ORC(Q_preheater, Q_boiler, W_turbine, Q_recuperator, Q_desuperheater, Q_condenser, ...
    W_pump_orc, W_pump_prod, dT_range_CT, dT_LMTD_preheater, dT_LMTD_boiler, dT_LMTD_recuperator, params)
%Heats and powers in Watts
%T in C

%Check heats & powers
if (Q_preheater < 0 || Q_boiler < 0)
    throw(MException('CapitalCost_SurfacePlant:NegativeBoilerHeat','Negative Boiler Heat'));
elseif (Q_recuperator < 0)
    throw(MException('CapitalCost_SurfacePlant:NegativeRecuperatorHeat','Negative Recuperator Heat')); 
elseif (W_turbine < 0)
    throw(MException('CapitalCost_SurfacePlant:NegativeTurbinePower','Negative Turbine Power')); 
elseif (Q_desuperheater > 0 || Q_condenser > 0)
    throw(MException('CapitalCost_SurfacePlant:PositiveCondenserHeat','Positive Condenser Heat')); 
elseif (W_pump_orc > 0)
    throw(MException('CapitalCost_SurfacePlant:PositiveOrcPumpPower','Positive ORC Pump Power')); 
elseif (W_pump_prod > 0)
    throw(MException('CapitalCost_SurfacePlant:PositiveProdPumpPower','Positive Prod Pump Power')); 
end

% Check temps
if (params.dT_approach < 0 || dT_range_CT < 0 || dT_LMTD_preheater < 0 || dT_LMTD_boiler < 0 || dT_LMTD_recuperator < 0)
    throw(MException('CapitalCost_SurfacePlant:NegativedT','Negative Temp Difference')); 
end

PPI_T_G = PPI('PPI_T-G', params.costYear); %2019 1.406
PPI_pump = PPI('PPI_Pump', params.costYear); %2019 1.617
PPI_HX = PPI('PPI_HX', params.costYear); %2019 1.797

% Calculate:
% C_T_G
%Regular fluid
S_T_fluid = 1.00; %Not CO2
C_T_G = 0.67 * PPI_T_G * (S_T_fluid*2830*(W_turbine/1e3)^0.745 + 3680*(W_turbine/1e3)^0.617);

% C_pump (ORC)
C_pump_orc_surface = 1750 * (1.34*-1*(W_pump_orc/1e3))^0.7;
S_pump_orc = 1.00; %Water
C_pump_orc = PPI_pump * S_pump_orc * C_pump_orc_surface;

% C_coolingTowers
TDC = 1;
[C_coolingTowers] = CapitalCost_CoolingTower(Q_desuperheater, Q_condenser, TDC, ...
    params.T_surface_air_C, params.dT_approach, dT_range_CT, params.costYear, params.coolingMode);

% C_heatExchanger
%dT_LMTD_HX
%U = 500/1000; %kW/m^2-K
%U = 500; %W/m^2-K
U = params.HX_overallHeatTransferCoefficient;
if (isnan(dT_LMTD_preheater) || dT_LMTD_preheater == 0)
    A_preheater = 0;
else
    A_preheater = Q_preheater / U / dT_LMTD_preheater;
end
A_boiler = Q_boiler / U / dT_LMTD_boiler;
A_HX = A_preheater + A_boiler;
C_heatExchanger = PPI_HX * (239*A_HX + 13400);

% C_recuperator
if (isnan(dT_LMTD_recuperator) || dT_LMTD_recuperator == 0 || Q_recuperator == 0)
    A_recuperator = 0;
    C_recuperator = 0;
else
    A_recuperator = Q_recuperator / U / dT_LMTD_recuperator;
    C_recuperator = PPI_HX * (239*A_recuperator + 13400);
end

% C_productionPump
C_pump_prod_lineshaft = 1750 * (1.34*-1*(W_pump_prod/1e3))^0.7 + 5750 * (1.34*-1*(W_pump_prod/1e3))^0.2;
S_pump_prod = 1.00; %Water
C_pump_prod = PPI_pump * S_pump_prod * C_pump_prod_lineshaft;

% THEN
% C_primaryEquipment
C_primaryEquipment = C_T_G + C_pump_orc + C_coolingTowers + C_heatExchanger + C_pump_prod + C_recuperator;

X_SE = 1.39;
X_CL = 0.58;
X_CM = 0.11;
X_ST = 0.00;
X_F = 0.04;
X_PC_sp = 1.15;
X_IC_sp = 1.12;

C_plant_TEC = X_SE * C_primaryEquipment;
C_plant_otherEquipment = C_plant_TEC - C_primaryEquipment;

C_plant_BEC = C_plant_TEC * (1 + X_CL + X_CM + X_ST + X_F);
C_plant_installation = C_plant_BEC - C_plant_TEC;

C_plant = X_PC_sp * X_IC_sp * C_plant_BEC;
C_plant_indirectContingency = C_plant - C_plant_BEC;

result.C_T_G = C_T_G;
result.C_pump_orc = C_pump_orc;
result.C_coolingTowers = C_coolingTowers;
result.C_heatExchanger = C_heatExchanger;
result.C_recuperator = C_recuperator;
result.C_pump_prod = C_pump_prod;
result.C_plant_otherEquipment = C_plant_otherEquipment;
result.C_plant_installation = C_plant_installation;
result.C_plant_indirectContingency = C_plant_indirectContingency;
result.C_plant = C_plant;

% calculate fractions
result.f_T_G = C_T_G / C_primaryEquipment;
result.f_pump_orc = C_pump_orc / C_primaryEquipment;
result.f_coolingTowers = C_coolingTowers / C_primaryEquipment;
result.f_heatExchanger = C_heatExchanger / C_primaryEquipment; 
result.f_pump_prod = C_pump_prod / C_primaryEquipment;

% specific capital cost of ORC
result.c_plant_noProdPump = C_plant / (W_turbine + W_pump_orc);
result.c_plant_inclProdPump = C_plant / (W_turbine + W_pump_orc + W_pump_prod);

end

