function maxBoilTemp_C = MaxSubcritORCBoilTemp(params)

% This function finds the saturated vapor pressure/temperature that has the
% largest entropy. This guarantees that with no superheating, the turbine
% will always only have vapor in it.

% Only need to do this once for a fluid
persistent fluid;
persistent temp;
if (strcmp(fluid,params.orcFluid))
    maxBoilTemp_C = temp;
    return;
end

% Different fluid
P_crit = CoolProp('PCRIT', "", 0, "", 0, params.orcFluid);

% assume min pressure is condensing at 0C
P_min = CoolProp('P', 'T', 0+273.15, 'Q', 1, params.orcFluid);
dP = (P_crit-P_min)/1000;

% find pressure with maximum entropy
vals = [];
for P = P_min:dP:P_crit
    s = CoolProp('SMASS', 'P', P, 'Q', 1, params.orcFluid);
    vals = [vals; P s];
end

max_s = max(vals(:,2));
row_index = find(vals(:,2) == max_s);
max_P = vals(row_index,1);
max_T = CoolProp('T', 'P', max_P, 'SMASS', max_s, params.orcFluid) - 273.15;

maxBoilTemp_C = max_T;

fluid = params.orcFluid;
temp = max_T;

end