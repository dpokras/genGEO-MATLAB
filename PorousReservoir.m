function result = PorousReservoir(P_f_initial, T_f_initial, m_dot, params)
   
% use transient characterization from Adams (2015)

% Assume that reservoir area, regardless of configuration, is 5-spot
A_reservoir = (params.res_length^2)/2;

%reservoirConfiguration = '5spot';
%reservoirConfiguration = 'Doublet';

time_seconds = 3600 * 24 * 365 * params.time_years; %s

% make sure temp gradient positive
if (params.dT_dz < 0)
    dT_dz = -params.dT_dz;
else
    dT_dz = params.dT_dz;
end
if (params.depth < 0)
    depth = -params.depth;
else
    depth = params.depth;
end

T_reservoir = depth*dT_dz + params.T_surface_rock_C;
P_reservoir = depth * 1000 * 9.81;

% evaluate viscosity and density
% VISCOSITY: Pa-s, DENSITY: kg/m^3
%from input properties
mu_inj_fluid = CoolProp('VISCOSITY', 'P', P_f_initial, 'T', T_f_initial+273.15, params.fluid);
rho_inj_fluid = CoolProp('DMASS', 'P', P_f_initial, 'T', T_f_initial+273.15, params.fluid);
cp_inj_fluid = CoolProp('CPMASS', 'P', P_f_initial, 'T', T_f_initial+273.15, params.fluid);
%from output properties
mu_prod_fluid = CoolProp('VISCOSITY', 'P', P_reservoir, 'T', T_reservoir+273.15, params.fluid);
rho_prod_fluid = CoolProp('DMASS', 'P', P_reservoir, 'T', T_reservoir+273.15, params.fluid);
cp_prod_fluid = CoolProp('CPMASS', 'P', P_reservoir, 'T', T_reservoir+273.15, params.fluid);
mu_fluid = (mu_inj_fluid + mu_prod_fluid)/2;
rho_fluid = (rho_inj_fluid + rho_prod_fluid)/2;
cp_fluid = (cp_inj_fluid + cp_prod_fluid)/2;

% reservoir impedance
if (params.res_length <= params.well_radius)
    RI = 0;
else
    if (contains(params.wellFieldType, '5spot'))
        A_c_rock = log( (4*params.res_length) / (2*params.well_radius*pi()) );
        RI = mu_fluid/rho_fluid/params.transmissivity * A_c_rock;
    elseif (strcmp(params.wellFieldType, 'Doublet') || strcmp(params.wellFieldType, 'Tungsten'))
        A_c_rock = log( (params.res_length) / (2*params.well_radius) );
        RI = mu_fluid/rho_fluid/params.transmissivity/pi() * A_c_rock;
    else
        throw(MException('PorousReservoir:unknownReservoirConfiguration','Unknown Reservoir Configuration'));
    end
end

% Calculate heat extracted
dT_initial = T_reservoir - T_f_initial;
res_energy = A_reservoir * params.thickness * params.rho_rock * params.C_rock * dT_initial;

% Model pressure transient (Figure 4.2, Adams (2015)), only for CO2 drying
% out
if (params.modelResPressureTransient == true && strcmp(params.fluid,'CO2'))
    R = params.res_length;
    nu_inj_fluid = mu_inj_fluid / rho_inj_fluid;
    tau = time_seconds * nu_inj_fluid / R^2;
    b = -0.4574 * (abs(R)/abs(params.depth))^-0.27;
    if (b < -0.9); b = -0.9; end
    if (b > -0.4); b = -0.4; end
    a = 1.0342 * exp(10.989*b);
    if (a < 0); a = 0; end
    if (a > 0.012); a = 0.012; end
    coeff = (a*tau^b + 1);
    RI = RI * coeff;
end

% Model temperature transient (Figure 4.5, Adams (2015))
% First-principles model for all fluids
% At time zero, output is the same as no temp depletion
if (params.modelResTemperatureDepletion == true && params.time_years > 0)
    % p1
    coeff1 = params.thickness * params.k_rock * rho_inj_fluid / m_dot / params.rho_rock / params.C_rock;
    p1 = -890.97*coeff1 - 1.30;
    % limit to regression
    if (p1 > -1.3); p1 = -1.3; end
    
    % p2
    p2 = 0.4095*exp(-0.7815*p1);
    
    % p3
    R = params.res_length;
    coeff2 = params.k_rock * R * rho_inj_fluid / params.rho_rock / cp_inj_fluid / m_dot;
    %p3 = 200.96*coeff2 - 0.144;
    % A more realistic, exponential relation is found by fitting the same data
    % with an exponential, instead of linear, curve.
    p3 = 1 - (1.4646*exp(-377.3*coeff2));
    % p3 can't be greater than 1 or less than 0.
    if (p3 < 0); p3 = 0; end
    if (p3 > 1); p3 = 1; end
    
    % psi
    % numerically integrate to get heat depletion
    %increments = 10;
    % have roughly one increment a year
    increments = round(params.time_years);
    if (increments < 1); increments = 1; end
    dt = time_seconds / increments;
    res_energy_extracted = 0;
    gamma = 1;
    for i = 1:increments
        res_energy_extracted = m_dot * cp_fluid * gamma * dT_initial * dt + res_energy_extracted;
        if (res_energy == 0)
            Psi = 0;
        else
            Psi = res_energy_extracted / res_energy;
        end
        gamma = DepletionCurve(Psi, p1, p2, p3);
    end
    
    % last gamma is res temp
else
    % gamma of 1 is undepleted reservoir
    gamma = 1;
    res_energy_extracted = m_dot * cp_fluid * gamma * dT_initial * time_seconds;
    if (res_energy <= 0)
        Psi = 0;
    else
        Psi = res_energy_extracted / res_energy;
    end
end

% Calculate final values
result.dP = m_dot * RI;
result.EndPressure = P_f_initial - result.dP;
result.EndTemp = gamma*dT_initial + T_f_initial;
result.EndEnthalpy = CoolProp('HMASS', 'P', result.EndPressure, 'T', result.EndTemp+273.15, params.fluid);
result.StartEnthalpy = CoolProp('HMASS', 'P', P_f_initial, 'T', T_f_initial+273.15, params.fluid);
result.Heat = m_dot * (result.EndEnthalpy - result.StartEnthalpy);

result.State = [params.res_length, result.dP/1e6, T_reservoir, T_reservoir];
result.Psi = Psi;
result.ReservoirImpedance = RI;
result.AvgFluidViscosity = mu_fluid;
result.AvgFluidDensity = rho_fluid;
result.AvgFluidSpecificHeat = cp_fluid;
   
end

