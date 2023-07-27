clear;
close all;

params = SimulationParameters;
well_radius = 0.2;
depth = 1000;
%fluid = 'CO2';
fluid = 'Water';

wellDiameter = 2*well_radius;
wellLength = depth;
successRate = 1;

[WellCost] = CapitalCost_Well(wellLength, wellDiameter, fluid, successRate, params);