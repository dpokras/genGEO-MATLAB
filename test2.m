dL = 1000
R = dL/(2*pi)
alpha = 30; %degrees
beta = linspace(0,2*pi,20) %radians
lambda = R - R*cos(beta)
z = lambda*sin(alpha*pi/180)

figure();
plot(beta*180/pi, lambda);
hold on;
plot(beta*180/pi, z);
% for i=1:6
% plot(r(i), z(i), ...
%     'Color', rand(1, 3), ...
%     'DisplayName',strcat('stream no. ',num2str(i)) ...
%     );
% end