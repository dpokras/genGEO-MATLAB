function [C_plant, f_T_G, f_pump_orc, f_coolingTowers, f_heatExchanger, f_pump_prod] = ...
    CapitalCost_SurfacePlant(Q_preheater, Q_boiler, W_turbine, Q_desuperheater, Q_condenser, ...
    W_pump_orc, W_pump_prod, T_ambient_C, dT_approach_CT, dT_range_CT, dT_LMTD_preheater, dT_LMTD_boiler)
%Heats and powers in kW
%T in C

%Check heats & powers
if (Q_preheater < 0 || Q_boiler < 0)
    throw(MException('CapitalCost_SurfacePlant:NegativeBoilerHeat','Negative Boiler Heat'));
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
if (dT_approach_CT < 0 || dT_range_CT < 0 || dT_LMTD_preheater < 0 || dT_LMTD_boiler < 0)
    throw(MException('CapitalCost_SurfacePlant:NegativedT','Negative Temp Difference')); 
end

CPI_T_G = 1.354; %2017
CPI_pump = 1.436; %2017
CPI_PE = 1.679; %2017
CPI_HX = 1.658; %2017

% Calculate:
% C_T_G
%Regular fluid
S_T_fluid = 1.00; %Not CO2
C_T_G = 0.67 * CPI_T_G * (S_T_fluid*2830*W_turbine^0.745 + 3680*W_turbine^0.617);

% C_pump (ORC)
C_pump_orc_surface = 1750 * (1.34*-1*W_pump_orc)^0.7;
S_pump_orc = 1.00; %Water
C_pump_orc = CPI_pump * S_pump_orc * C_pump_orc_surface;

% C_coolingTowers
a_cool = 5.58e3;
b_cool = 0;
c_cool = -1.77e1;
d_cool = 1.96e2;
c_cooling_wet = a_cool*(1/dT_approach_CT) + b_cool*(T_ambient_C+273.15) + c_cool*(T_ambient_C+273.15)/dT_approach_CT + d_cool*(1/(dT_approach_CT + dT_range_CT)); 
a_cond = 4.08e3;
b_cond = -1.54e-2;
c_cond = -1.24e1;
d_cond = 0;
c_condensing_wet = a_cond*(1/dT_approach_CT) + b_cond*(T_ambient_C+273.15) + c_cond*(T_ambient_C+273.15)/dT_approach_CT + d_cond*(1/(dT_approach_CT + dT_range_CT)); 
B_cooling = 1.00; %Closed-circuit R245fa cooler
B_condensing = 1.00; %Closed-circuit R245fa condenser
C_BAC = CPI_PE * (B_cooling*c_cooling_wet*-1*Q_desuperheater + B_condensing*c_condensing_wet*-1*Q_condenser);
C_coolingTowers = C_BAC * (-1*(Q_desuperheater + Q_condenser)/1000);

% C_heatExchanger
%dT_LMTD_HX
U = 750/1000; %kW/m^2-K
A_preheater = Q_preheater / U / dT_LMTD_preheater;
A_boiler = Q_boiler / U / dT_LMTD_boiler;
A_HX = A_preheater + A_boiler;
C_heatExchanger = CPI_HX * (239*A_HX + 13400);

% C_productionPump
C_pump_prod_lineshaft = 1750 * (1.34*-1*W_pump_prod)^0.7 + 5750 * (1.34*-1*W_pump_prod)^0.2;
S_pump_prod = 1.00; %Water
C_pump_prod = CPI_pump * S_pump_prod * C_pump_prod_lineshaft;

% THEN
% C_primaryEquipment
C_primaryEquipment = C_T_G + C_pump_orc + C_coolingTowers + C_heatExchanger + C_pump_prod;

X_SE = 1.39;
X_CL = 0.58;
X_CM = 0.11;
X_ST = 0.00;
X_F = 0.04;
X_PC_sp = 1.15;
X_JC_sp = 1.10;
X_IC_sp = 1.12;

C_plant_TEC = X_SE * C_primaryEquipment;

C_plant_BEC = C_plant_TEC * (1 + X_CL + X_CM + X_ST + X_F);

C_plant = X_PC_sp * X_JC_sp * X_IC_sp * C_plant_BEC;

% calculate fractions
f_T_G = C_T_G / C_primaryEquipment;
f_pump_orc = C_pump_orc / C_primaryEquipment;
f_coolingTowers = C_coolingTowers / C_primaryEquipment;
f_heatExchanger = C_heatExchanger / C_primaryEquipment; 
f_pump_prod = C_pump_prod / C_primaryEquipment;

end

