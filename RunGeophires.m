function result = RunGeophires(reservoirImpedance, reservoirFluidDensity, m_dot_well, numInjWells, numProdWells, params, pathToScripts)

%Check that flowrate is between 1 and 500
if (m_dot_well < 1 || m_dot_well > 500)
    throw(MException('RunGeophires:Bad_m_dot','Bad Mdot'));
end

%impedance is GPa-s/m^3, but we use Pa-s/kg
impedance = reservoirImpedance * reservoirFluidDensity / 1e9;

if (impedance < 1e-4 || impedance > 1e4)
    throw(MException('RunGeophires:Bad_impedance',strcat(['Bad Impedance of ' num2str(impedance) ', must be between 1e-4 and 1e4'])));
end

gradient = params.dT_dz * 1000;
depth = params.depth / 1000;
wellDiameter = params.well_radius*2/0.3048*12;
rockDensity = params.rho_rock;
rockHeatCapacity = params.C_rock;
rockThermalConductivity = params.k_rock;
eta_pump = params.eta_pump;
capacityFactor = params.CapacityFactor;
surfaceRockTemp = params.T_surface_rock_C;
surfaceAirTemp = params.T_surface_air_C;
plantLifetime = params.Lifetime;
discountRate = params.discountRate;

%impedance = 0.1; %GPa-s/m^3
%m_dot = 10; %kg/s
%gradient = 35; %C/km
%depth = 2.5; %km
%wellDiameter = 12.25; %inches
%rockDensity = 2650; %kg/m^3
%rockHeatCapacity = 1000; %J/kg/K
%rockThermalConductivity = 2.1; %W/m^2-K
%eta_pump = 0.9;
%capacityFactor = 0.85;
%surfaceRockTemp = 15; %C
%surfaceAirTemp = 15; %C
%plantLifetime = 25; %years
%discountRate = 0.096;

%pathToScripts = 'C:\Users\Shackleton\Desktop\Dropbox\PhD\GEOPHIRES\';

data = fileread(fullfile(pathToScripts, 'input_generic.txt'));

data = strrep(data, '@impedance', num2str(impedance));
data = strrep(data, '@prodwell', num2str(m_dot_well));
data = strrep(data, '@numInjWells', num2str(numInjWells));
data = strrep(data, '@numProdWells', num2str(numProdWells));
data = strrep(data, '@gradient', num2str(gradient));
data = strrep(data, '@depth', num2str(depth));
data = strrep(data, '@wellDiameter', num2str(wellDiameter));
data = strrep(data, '@rockDensity', num2str(rockDensity));
data = strrep(data, '@rockHeatCapacity', num2str(rockHeatCapacity));
data = strrep(data, '@rockThermalConductivity', num2str(rockThermalConductivity));
data = strrep(data, '@circulationPumpEfficiency', num2str(eta_pump));
data = strrep(data, '@capacityFactor', num2str(capacityFactor));
data = strrep(data, '@surfaceRockTemp', num2str(surfaceRockTemp));
data = strrep(data, '@surfaceAirTemp', num2str(surfaceAirTemp));
data = strrep(data, '@plantLifetime', num2str(plantLifetime));
data = strrep(data, '@discountRate', num2str(discountRate));

fid = fopen(fullfile(pathToScripts, 'input_specific.txt'), 'wt');
fprintf(fid,'%s',data);
fclose(fid);

[status,cmdout] = system(strcat([ 'python ' fullfile(pathToScripts, 'GEOPHIRESv2.py') ]));

% Get Pump dP from result
data = fileread(fullfile(pathToScripts, 'HDR.txt'));
data = splitlines(data);
for i = 1:size(data,1)
    if (contains(data(i),'Average total geofluid pressure drop (kPa)'))
        break;
    end
end
output = split(strtrim(data(i)));
% fifth value is the LCOE $/kWh
dP_pump_kPa = str2num(cell2mat(output(7)));


% Get LCOE from result in cent/kWh
LCOE_Cent_kWh = getValueFromLine('Electricity breakeven price', fullfile(pathToScripts, 'HDR.txt'));
result.LCOE = LCOE_Cent_kWh / 100 / 1000; %to get $/Wh

% CAPEX in M$/yr
CapitalCost = getValueFromLine('Total capital costs', fullfile(pathToScripts, 'HDR.txt'));
result.CapitalCost = CapitalCost * 1e6; % to get $

% Drilling Costs (M$/yr)
DrillingCost = getValueFromLine('Drilling and completion costs', fullfile(pathToScripts, 'HDR.txt'));
result.C_wells = DrillingCost * 1e6;

% Surface Plant Costs (M$/yr)
SurfaceCost = getValueFromLine('Surface power plant costs', fullfile(pathToScripts, 'HDR.txt'));
result.C_surfacePlant = SurfaceCost * 1e6;

% Field Gathering Cost (M$/yr)
FieldCost = getValueFromLine('Field gathering system costs', fullfile(pathToScripts, 'HDR.txt'));
result.C_wellfield = FieldCost * 1e6;

% Exploration Cost (M$/yr)
ExplorationCost = getValueFromLine('Exploration costs', fullfile(pathToScripts, 'HDR.txt'));
result.C_exploration = ExplorationCost * 1e6;

% O&M Cost (M$/yr)
OMCost = getValueFromLine('Total operating and maintenance costs', fullfile(pathToScripts, 'HDR.txt'));
result.C_OM = OMCost * 1e6;


% Get power from result
data = fileread(fullfile(pathToScripts, 'HDR.txt'));
data = splitlines(data);
for i = 1:size(data,1)
    if (contains(data(i),'YEAR') && contains(data(i),'THERMAL'))
        break;
    end
end
% Now we are at the output yearly results
% skip three lines plus header
i = i + 4;
output = split(strtrim(data(i)));
%3rd - Production Temp
T_prod_surface_C = str2double(cell2mat(output(3)));
%4th - Pump Power
W_pump_prod = str2double(cell2mat(output(4)));
%5th - Net Power
W_net = str2double(cell2mat(output(5)));

result.T_prod_surface_C = T_prod_surface_C;
result.dP_pump = dP_pump_kPa * 1e3; % Pa
result.W_pump_prod = W_pump_prod * 1e6; %Watts
result.W_net = W_net * 1e6; % to get Watts
end

function value = getValueFromLine(str, path)
    data = fileread(path);
    data = splitlines(data);
    for i = 1:size(data,1)
        if (contains(data(i),str))
            break;
        end
    end
    
    match = regexp(data(i), '[.0-9]', 'match');
    matchstr = strjoin(match{1},'');
    value = str2double(matchstr);
    %output = split(strtrim(data(i)));
    % fifth value is the LCOE $/kWh
    %value = str2double(cell2mat(output(5)));
end
