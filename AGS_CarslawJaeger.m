% ------------------------------------------------------------------------
% ------------ Solution to Carlow & Jaeger (1959) Figure 41 --------------
% ------------------------------------------------------------------------

clear;
close all;
clc;

% Specify Temperature drawdown Tolerance
tolerance = 0.05;
% Figure 1 Plot for Tabled Values (Jaeger 1956)---------------------------
R = xlsread('data\Extracted Data Carslow.xlsx','B6:U209');
td = xlsread('data\Extracted Data Carslow.xlsx','W4:AQ4');
td = td';
td(any(isnan(td), 2)) = [];
legendCell(1) = arrayfun(@(x) sprintf('t = %.1f minutes',x),td(1),'un',0)
legendCell(2) = arrayfun(@(x) sprintf('t = %.1f   hours',x),td(2),'un',0)
legendCell(3) = arrayfun(@(x) sprintf('t = %.1f   hours',x),td(3),'un',0)
legendCell(4) = arrayfun(@(x) sprintf('t = %.1f   days',x),td(4),'un',0)
legendCell(5) = arrayfun(@(x) sprintf('t = %.1f   days',x),td(5),'un',0)
legendCell(6) = arrayfun(@(x) sprintf('t = %.1f   days',x),td(6),'un',0)
legendCell(7) = arrayfun(@(x) sprintf('t = %.1f days',x),td(7),'un',0)
legendCell(8) = arrayfun(@(x) sprintf('t = %.2f years',x),td(8),'un',0)
legendCell(9) = arrayfun(@(x) sprintf('t = %.2f years',x),td(9),'un',0)
legendCell(10) = arrayfun(@(x) sprintf('t = %.1f   years',x),td(10),'un',0)
legendCell(11) = arrayfun(@(x) sprintf('t = %.1f years',x),td(11),'un',0)
% Create vector showing temperature drawdown tolerance
ylimit = 100;
X = linspace(0.1,ylimit,100);
v_V = zeros(1,100);
for j = 1:100
    v_V(j) = tolerance;
end
color_wheel = ["#036382","#027ea6","#09aee3","#40b7dd","#1479d9","#1f90be","#6BA7CC","#AEDBF0","#CBF1FA","#E2FCFF"]
j=1;
for i = 1:2:length(td)*2-2
    RR = R(:,i:i+1);
    RR(any(isnan(RR), 2), :) = [];
    t = RR(:,2); % Retrieve temperature values
    r = RR(:,1); % Retrieve radial distance values
    tt = smooth(t);
    rr = smooth(r);
    figure(1)
    semilogx(rr,tt,'color',color_wheel(j),'linewidth',1.4);
    j = j + 1;
    grid on
    xlim([min(r),ylimit]) % Specifying limits on x-axis
    xlabel('Radial Distance from Center of the Well, r [m]', 'FontSize',11) % x-axis label
    ylabel('\Gamma_{r} [ - ]', 'FontSize',11) % y-axis label
    set(gca,'FontSize',11)
    hold on
end
figure(1)
extrap_r = 38.348564;
year_30_extrap = plot(extrap_r,0.05,'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black');
tolerance = semilogx(X,v_V,'--k');
lgd = legend(legendCell,'numcolumns',2);
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 11;
lgd.FontSize = 11;

saveas(figure(1),'images\Carslaw_Jaeger_Fig41_1959.fig');