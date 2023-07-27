% x = linspace(-2*pi,2*pi,10);
% y = linspace(0,4*pi,10);
% [X,Y] = meshgrid(x,y);
% Z = sin(X)+cos(Y);
% contour(X,Y,Z)
for l = 5000:2000:15000
    filename = strcat('dpo_opt_config1_res',num2str(l),'_05_06_23.csv')
    M = csvread(filename);
end