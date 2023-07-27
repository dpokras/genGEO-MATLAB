function ff = FrictionFactor(wellRadius, P, h, m_dot, params)
rho_fluid = CoolProp.PropsSI('DMASS', 'P', P, 'HMASS', h, params.fluid);
mu = CoolProp.PropsSI('VISCOSITY', 'P', P, 'HMASS', h, params.fluid);

A_c = pi * wellRadius.^2;
V = m_dot / A_c ./ rho_fluid;

% assert(V(1)<4,'velocity is too high. (V > 4m/s)')

% Relative Roughness
rr = params.epsilon/(2*wellRadius);

% Reynolds Number
Re = rho_fluid.*V*(2*wellRadius)./mu;

% Use Haaland (1983) Approximation
c1 = -1.8*log10((6.9./Re) + (rr/3.7)^1.11);
ff = (1./c1).^2;

end