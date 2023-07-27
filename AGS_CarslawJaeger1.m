% ------------------------------------------------------------------------
% ------------ Solution to Carlow & Jaeger (1959) Figure 41 --------------
% ------------------------------------------------------------------------

format short
clear;
close all;

% Specify Temperature drawdown Tolerance
tolerance = 0.05;
% Figure 1 Plot for Tabled Values (Jaeger 1956)---------------------------
R = xlsread('data\Extracted Data Carslow.xlsx','B6:U209');
td = xlsread('data\Extracted Data Carslow.xlsx','B4:U4');
td = td';
td(any(isnan(td), 2)) = [];
legendCell = arrayfun(@(x) sprintf('t = %.1d years',x),td,'un',0);

% Create vector showing temperature drawdown tolerance
X = linspace(1,max(R(:,19)),100);
v_V = zeros(1,100);
for j = 1:100
    v_V(j) = tolerance;
end


RR = R(:,1:2);
RR(any(isnan(RR), 2), :) = [];
t = RR(:,2); % Retrieve temperature values
r = RR(:,1); % Retrieve radial distance values
tt = smooth(t);
rr = smooth(r);
figure(1)
semilogx(rr,tt,'-k');
grid on
xlim([min(r),max(r)]) % Specifying limits on x-axis
xlabel('r/a [ - ]','interpreter','latex'); % x-axis label
ylabel('\Gamma_{r}');    % y-axis label

hold on
semilogx(r,t,'r')

figure(1)
tolerance = semilogx(X,v_V,'--r');



saveas(figure(1),'images\Carslaw_Jaeger_1.fig');