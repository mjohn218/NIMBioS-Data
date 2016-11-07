clear
clc

load('minDt.mat');
load('time.mat');
load('distance.mat');

figure(2)
surf(time,distance,minDt');
axis([0 100 0 6]);
set(gca, 'FontSize',16);
xlabel('time(s)', 'FontSize',16);
ylabel('distance along major axis({\mu}m)','FontSize',16);
title('Density of minDt on Membrane','FontSize',18);
c=colorbar;
ylabel(c,'N_{minDt}(t)');