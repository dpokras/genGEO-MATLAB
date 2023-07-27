function result = semi_analytic_well(P_f_initial, T_f_initial, T_e_initial, dz_total, dr_total, m_dot, params)
%% Begin Iterations
N_dz = 25; % segmentation of the well
m = params.n_streams;
N_dt = size(params.time_years,2);

%% Calculate time dimension
t = permute(repmat(params.time_years.', [1, m, N_dz]), [3, 2, 1]);
time_seconds = 3600 * 24 * 365 * t; %s
g=9.81; %m/s^2


%% Calculate gradual decreasing inner radius for vertical wells only.
if dr_total == 0 & dz_total < 1 % injection well; calculating vertical well radii
    wellRadius = ones(N_dz,m,N_dt).* linspace(params.well_radius, params.side_stream_radius, N_dz).';
elseif dr_total == 0 & dz_total > 1 %production well; calculating vertical well radii
    wellRadius = ones(N_dz,m,N_dt).* flip(linspace(params.well_radius, params.side_stream_radius, N_dz)).';
else % reservoir; calculating lateral well radius
    wellRadius = ones(N_dz,m,N_dt) .* params.side_stream_radius;
end

%% dz_total and dr_total are arrays

dz = dz_total/N_dz; %m
dr = dr_total/N_dz; %m

dL = (dz.^2 + dr.^2).^0.5; %m
A_c = pi .* wellRadius.^2; %m^2
P_c = pi .* 2 .* wellRadius; %m

z = zeros(N_dz,m,N_dt);
r = zeros(N_dz,m,N_dt);
L = zeros(N_dz,m,N_dt);
T_f = zeros(N_dz,m,N_dt);
T_e = zeros(N_dz,m,N_dt);
P = zeros(N_dz,m,N_dt);
h = zeros(N_dz,m, N_dt);
rho_fluid = zeros(N_dz,m,N_dt);
q = zeros(N_dz,m,N_dt);
Cp_fluid = zeros(N_dz,m,N_dt);
v = zeros(N_dz,m,N_dt);
delta_P_loss = zeros(N_dz,m,N_dt);

T_f(1,:,:) = ones(1, m, N_dt) .* T_f_initial; %C
T_e(1,:,:) = ones(1, m, N_dt) .* T_e_initial; %C
P(1,:,:) = ones(1, m, N_dt) .* P_f_initial; %z(1) * 10000; %Pa


h(1,:,:) = CoolProp.PropsSI('HMASS', 'P', P(1,:,:), 'T', T_f(1,:,:)+273.15, params.fluid);
rho_fluid(1,:,:) = CoolProp.PropsSI('DMASS', 'P', P(1,:,:), 'T', T_f(1,:,:)+273.15, params.fluid);
v(1,:) = m_dot./A_c(1)./rho_fluid(1,:);

% Use Colebrook-white equation for wellbore friction loss.
% Calculate frictionFactor for first element and assume constant in remainder
ff = FrictionFactor(wellRadius(1), P(1,:,:), h(1,:,:), m_dot, params);


alpha_rock = params.k_rock / params.rho_rock / params.C_rock; %dim

% HeatLoss
% For heatloss calculation
t_d = alpha_rock.*time_seconds./(wellRadius).^2; %dim
if (t_d < 2.8)
    beta = ((pi*t_d)^-0.5 + 0.5 - 0.25*(t_d/pi)^0.5 + 0.125*t_d);
else
    beta = (2./(log(4.*t_d)-2.*0.58) - 2.*0.58./(log(4.*t_d)-2.*0.58).^2);
end

for i=2:N_dz+1
    
    L(i,:,:)=L(i-1,:,:)+dL;

    if dz_total == 0
        R = dr_total/(2*pi);
        z(i,:,:) = -(R-R.*cos(L(i,:,:)./R)).*sin(params.angle*pi/180);
    else
        z(i,:,:)=z(i-1,:,:)+dz;
    end

    r(i,:,:)=r(i-1,:,:)+dr;
    
    % far-field rock temp
    T_e(i,:,:) = T_e_initial - z(i,:,:) * params.dT_dz;
    % fluid velocity
    v(i,:,:) = m_dot/A_c(i-1)./rho_fluid(i-1,:,:); %m/s
    
    % Calculate Pressure
    delta_P_loss(i,:,:) = ff .* dL ./ (2.*wellRadius(i-1)) .* rho_fluid(i-1,:,:) .* v(i,:,:).^2 ./ 2; %Pa
    rho_fluid(i,:,:) = CoolProp.PropsSI('DMASS', 'P', P(i-1,:,:), 'HMASS', h(i-1,:,:), params.fluid);
    adj_var = 3e3;
    if any(isinf(rho_fluid(i,:,:))) == 1
        rho_fluid(i,:,:) = CoolProp.PropsSI('DMASS', 'P', P(i-1,:,:)+adj_var, 'HMASS', h(i-1,:,:), params.fluid);
    end
    P(i,:,:) = P(i-1,:,:) - rho_fluid(i,:,:).*g.*dz - delta_P_loss(i,:,:);

    % Throw exception if below saturation pressure of water at previous temperature
    if (strcmp(params.fluid,'Water'))
        P_sat = CoolProp.PropsSI('P', 'T', T_f(i-1,:,:) + 273.15, 'Q', 0, params.fluid);
        if (P(i,:) < P_sat)
            throw(MException('SemiAnalytic:BelowSaturationPressure',strcat(['Below saturation pressure of water at ' num2str(z(i,:)) 'm !'])));
        end
    end
    
    h_noHX = h(i-1,:,:) - g.*dz;
    T_noHX = CoolProp.PropsSI('T', 'P', P(i,:,:), 'HMASS', h_noHX, params.fluid) - 273.15;

    Cp_fluid(i,:,:) = CoolProp.PropsSI('CPMASS', 'P', P(i,:,:), 'HMASS', h_noHX, params.fluid);
    % sometimes a random number might break in the CoolProp package.
    % Usually close to the critical point. this is a problem with the
    % package and not a physical problem. The pressure was simply moved by
    % 2 kPa
    if any(isinf(T_noHX)) == 1
        T_noHX = CoolProp.PropsSI('T', 'P', P(i,:,:) + adj_var, 'HMASS', h_noHX, params.fluid) - 273.15;
        Cp_fluid(i,:,:) = CoolProp.PropsSI('CPMASS', 'P', P(i,:,:) + adj_var, 'HMASS', h_noHX, params.fluid);
    end
    
    %Find Fluid Temp
    if (params.useWellboreHeatLoss == false)
        T_f(i,:,:) = T_noHX;
        h(i,:,:) = h_noHX;
        q(i,:,:) = 0;
    else
        % See Zhang, Pan, Pruess, Finsterle (2011). A time-convolution
        % approach for modeling heat exchange between a wellbore and
        % surrounding formation. Geothermics 40, 261-266.
        x = dL.*P_c(i-1,:,:).*params.k_rock.*beta(i-1,:,:)./wellRadius(i-1,:,:);
        y = m_dot.*Cp_fluid(i,:,:);
        %ratio = y/x;
        if (isinf(x))
            T_f(i,:,:) = T_e(i,:,:);
        else
            T_f(i,:,:) = (y.*T_noHX + x.*T_e(i,:,:)) ./ (x + y);
        end
        q(i,:,:) = y.*(T_noHX - T_f(i,:,:));
        h(i,:,:) = CoolProp.PropsSI('HMASS', 'P', P(i,:,:), 'T', T_f(i,:,:)+273.15, params.fluid);
    end
end

result.State = [z(end,:,:), P(end,:,:)/1e6, T_f(end,:,:), T_e(end,:,:)];

result.stream_depth = z(end,:,end);
result.EndTemp = T_f(end,:,:);
result.EndPressure = P(end,:,:);
result.EndEnthalpy = h(end,:,:);
result.EndDensity = rho_fluid(end,:,:);
result.Heat = -1*sum(q(end,:,:));
result.dP_loss = sum(delta_P_loss(end,:,:));
result.dL = dL;
result.N_dz = N_dz;
result.stream_depth = z(:,:,end);
result.res_length = r(:,:,end);
result.length = L(:,:,end);
result.Temp = T_f;
result.Pressure = P;
result.Enthalpy = h;
result.Density = rho_fluid;
result.wellRadius = wellRadius(:,end,end);
end

