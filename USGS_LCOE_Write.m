clear all;
close all;

Xtype = 'Transmissivity';

%xlsxFile = 'data\Porous_Feb2020.xlsx';
%xlsxFile = 'data\Porous_Jun2020.xlsx';
%xlsxFile = 'data\Porous_Jul2020.xlsx';
%xlsxFile = 'data\Porous_Aug2020.xlsx';
xlsxFile = 'data\Porous_Oct2020.xlsx';
%xlsxFile = 'data\Conduction_Jun2020.xlsx';
co2Sheet = 'CO2';
waterSheet = 'Water (R245fa)';
%waterSheet = 'Water (R600a)';

%Get LCOE
[X_LCOE_Water,Y_LCOE_Water,LCOE_Water,LCOE_Combined_Mask_Water] = GetContourData('Water', false, Xtype, 'LCOE_Greenfield', xlsxFile, co2Sheet, waterSheet);
[X_LCOE_CO2,Y_LCOE_CO2,LCOE_CO2,LCOE_Combined_Mask_CO2] = GetContourData('CO2', false, Xtype, 'LCOE_Brownfield', xlsxFile, co2Sheet, waterSheet);

%Get Power
[X_Power_Water,Y_Power_Water,Power_Water,Power_Combined_Mask_Water] = GetContourData('Water', false, Xtype, 'Power', xlsxFile, co2Sheet, waterSheet);
[X_Power_CO2,Y_Power_CO2,Power_CO2,Power_Combined_Mask_CO2] = GetContourData('CO2', false, Xtype, 'Power', xlsxFile, co2Sheet, waterSheet);

[usgs_data, usgs_headers, usgs_complete] = xlsread('data\USGS_Table_1.xlsx', 'Table_1');
%area is column 3 (1 without row headers) (1 acre = 0.004047 km^2)
%area = usgs_data(:,1)*0.004047;
%thickness is column 4 (2 without row headers)
thickness = usgs_data(:,2) * 0.3048;
%porosity is column 5 (really 3)
porosity = usgs_data(:,3);
%perm is column 7 (really 5)
perm = usgs_data(:,5);
transmissivity = thickness .* perm;
%depth is column 9 (really 7)
depth = usgs_data(:,7) * 0.3048;
%co2 density is column 8 (really 6)
co2density = usgs_data(:,6);
%technical co2 storage of formation is column 11 (really 9)
techco2storage = usgs_data(:,9) * 1000 * 1000 * 1000; %to get kg

%calculate co2 for 1 km^2
efficiencyfactor = 0.36;
volume = 1000 * 1000 * thickness .* porosity * efficiencyfactor;
fivespotmass = volume .* co2density;
% area is now number of 5spots worth of co2 that can be used
area = techco2storage ./ fivespotmass;

% calculate LCOE
% LCOE Factor (LCOE_LAZARD -> 1, LCOE_LIKELY -> 0.701
lcoe_factor = 0.701;
lcoe_CO2 = lcoe_factor * interp2(X_LCOE_CO2,Y_LCOE_CO2,LCOE_CO2,log10(transmissivity),depth);
lcoe_Water = lcoe_factor * interp2(X_LCOE_Water,Y_LCOE_Water,LCOE_Water,log10(transmissivity),depth);

% calculate Power
pow_CO2 = interp2(X_Power_CO2,Y_Power_CO2,Power_CO2,log10(transmissivity),depth);
pow_Water = interp2(X_Power_Water,Y_Power_Water,Power_Water,log10(transmissivity),depth);

% calculate total power
total_power_CO2 = pow_CO2 .* area;
total_power_Water = pow_Water .* area;

% extract only first two columns
newheader = usgs_headers(:,1:2);
newlcoe1 = cat(1,'LCOE_CO2 ($/MWh)', 'LCOE_CO2 ($/MWh)', num2cell(lcoe_CO2));
newlcoe2 = cat(1,'LCOE_Water ($/MWh)', 'LCOE_Water ($/MWh)',num2cell(lcoe_Water));
newpow1 = cat(1,'Power_CO2 (MWe)', 'Power_CO2 (MWe)', num2cell(pow_CO2));
newpow2 = cat(1,'Power_Water (MWe)', 'Power_Water (MWe)',num2cell(pow_Water));
newtotpow1 = cat(1,'Total_Power_CO2 (MWe)', 'Total_Power_CO2 (MWe)', num2cell(total_power_CO2));
newtotpow2 = cat(1,'Total_Power_Water (MWe)', 'Total_Power_Water (MWe)',num2cell(total_power_Water));

newdata = [usgs_complete newlcoe1 newlcoe2 newpow1 newpow2 newtotpow1 newtotpow2];
orig_col_number = size(usgs_complete,2);

%remove redundant header row
newdata(1,:) = [];

%sort rows
data_sorted_lcoeco2 = sortrows(newdata(2:end,:),orig_col_number + 1);
data_sorted_lcoewater = sortrows(newdata(2:end,:),orig_col_number + 2);

%add header row back
data_sorted_lcoeco2 = [newdata(1,:); data_sorted_lcoeco2];
data_sorted_lcoewater = [newdata(1,:); data_sorted_lcoewater];

%add power cumulatively (CO2)
clear cumul_pow_co2;
cumul_pow_co2(1) = 0;
for i=2:size(data_sorted_lcoeco2,1)
    cumul_pow_co2(i) = cumul_pow_co2(i-1) + cell2mat(data_sorted_lcoeco2(i,orig_col_number + 5));
end
cumul_pow_co2 = num2cell(transpose(cumul_pow_co2));
cumul_pow_co2(1) = {'Cumulative Power CO2 (MWe)'};
data_sorted_lcoeco2 = [data_sorted_lcoeco2 cumul_pow_co2];
%add power cumulatively (water)
clear cumul_pow_water;
cumul_pow_water(1) = 0;
for i=2:size(data_sorted_lcoewater,1)
    cumul_pow_water(i) = cumul_pow_water(i-1) + cell2mat(data_sorted_lcoewater(i,orig_col_number + 6));
end
cumul_pow_water = num2cell(transpose(cumul_pow_water));
cumul_pow_water(1) = {'Cumulative Power Water (MWe)'};
data_sorted_lcoewater = [data_sorted_lcoewater cumul_pow_water];

%add co2 cumulatively (CO2)
clear cumul_co2;
cumul_co2(1) = 0;
for i=2:size(data_sorted_lcoeco2,1)
    cumul_co2(i) = cumul_co2(i-1) + cell2mat(data_sorted_lcoeco2(i,11));
end
cumul_co2 = num2cell(transpose(cumul_co2));
cumul_co2(1) = {'Cumulative CO2 Needed (Mt)'};
data_sorted_lcoeco2 = [data_sorted_lcoeco2 cumul_co2];

%plot capacity
x_co2 = cell2mat(data_sorted_lcoeco2(2:end,18))/1000;
y_co2 = cell2mat(data_sorted_lcoeco2(2:end,12));
y_co2_needed = cell2mat(data_sorted_lcoeco2(2:end,19))/1000;
x_water = cell2mat(data_sorted_lcoewater(2:end,18))/1000;
y_water = cell2mat(data_sorted_lcoewater(2:end,13));
hold on;
stairs(x_co2,y_co2,'Color','#F79646','LineWidth',2);
stairs(x_water,y_water,'Color','k','LineWidth',2);
yyaxis right;
plot(x_co2,y_co2_needed);
ylabel(gca,'CO_2 Required [Gt]');
set(gca, 'Ylim', [0 500]);
yyaxis left;
stem(x_co2,y_co2,'LineStyle',':','Marker','none','Color','#F79646','LineWidth',1.25);
stem(x_water,y_water,'LineStyle',':','Marker','none','Color','k','LineWidth',1.25);
hold off
ylim([0 300]);
xlim([0 120]);
set(gcf, 'Position', [100, 100, 1000, 600]);
set(gca, 'FontSize', 16);
xlabel('Electric Generating Capacity [GWe]');
ylabel('LCOE_{LIKELY} [2019$/MW-h]');
legend({'CO_2','Water','CO_2 Required'});
grid on;

saveas(gcf,'images\USGS_Capacity.png');

%write to file
delete 'data\USGS_Table_MatlabResults.xlsx';
writecell(data_sorted_lcoeco2,'data\USGS_Table_MatlabResults.xlsx', 'Sheet', 'Matlab_Output_sortedco2');
writecell(data_sorted_lcoewater,'data\USGS_Table_MatlabResults.xlsx', 'Sheet', 'Matlab_Output_sortedwater');
