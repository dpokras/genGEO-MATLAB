clear all;
close all;

discountRate = 0:1:12;
OMRate = 0:1:10;
LT = 25;

[X,Y] = meshgrid(discountRate,OMRate);

for i = 1:size(discountRate,2)
    d = discountRate(i);
    if (d == 0)
        CRF = 1/LT;
    else
        CRF = (d/100)*(1+(d/100))^LT / ((1+(d/100))^LT-1);
    end
    for j = 1:size(OMRate,2)
        o = OMRate(j);
        CF = 0.85;
        Z85(j,i) = (CRF + (o/100)) / CF / 8760 * 1e6;
        CF = 0.9;
        Z90(j,i) = (CRF + (o/100)) / CF / 8760 * 1e6;
        CF = 0.95;
        Z95(j,i) = (CRF + (o/100)) / CF / 8760 * 1e6;
    end
end


hold on;
[M85,c85] = contour(X,Y,Z85,'--','Color','#888','DisplayName','CF=0.85');
c85.LevelList = [10, 20];
c85.LabelSpacing = 200;
c85.ShowText = 1;
%clabel(M85,c85,'manual','Color',[0.5 0.5 0.5]);
clabel(M85,c85,'Color',[0.5 0.5 0.5]);
[M90,c90] = contour(X,Y,Z90,'-','Color','b','DisplayName','CF=0.90');
c90.LevelList = 2:2:30;
c90.LabelSpacing = 200;
c90.ShowText = 1;
%clabel(M90,c90,'manual');
[M95,c95] = contour(X,Y,Z95,':','Color','#888','DisplayName','CF=0.95');
c95.LevelList = [10, 20];
c95.LabelSpacing = 200;
c95.ShowText = 1;
%clabel(M95,c95,'manual','Color',[0.5 0.5 0.5]);
clabel(M95,c95,'Color',[0.5 0.5 0.5]);

scatter(9.6,4.5,'o','DisplayName','Lazard','LineWidth',1.5);
scatter(4,5.5,'^','DisplayName','Likely (genGEO)','LineWidth',1.5);
scatter(0,1.5,60,'h','DisplayName','O&M-Only (genGEO)','LineWidth',1.5);
scatter(7,3,'d','DisplayName','GETEM','LineWidth',1.5);
scatter(7,6,'s','DisplayName','GEOPHIRES','LineWidth',1.5);
hold off;
xlabel('Weighted Average Cost of Capital, WACC (Discount Rate) [%]');
xticks(0:1:12);
ylabel('O&M Percent of Capital Cost, F_{O&M}  [%/yr]');
title('$$RF = \frac{LCOE}{SpCC}$$  ( x10$^{-6}$) [1/hr]','Interpreter','latex');
legend();
grid on;
box on;

set(gcf,'Position',[50,50,600,450]);

fileName = strcat(['images/lcoe_coefficient.png']);
saveas(gcf,fileName);