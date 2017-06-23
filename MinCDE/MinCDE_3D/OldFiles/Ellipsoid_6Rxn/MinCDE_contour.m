%% Reinitialize
clear
clc
clf

%% Import
load('minDt.mat');
load('EminDt.mat');
load('time.mat');
load('distance.mat');

%% Plot
% MinDt Plot
figure(2)
surf(time,distance,minDt');
axis([0 1200 0 6]);
set(gca, 'FontSize',12);
xlabel('time(s)', 'FontSize',12);
ylabel('distance along major axis({\mu}m)','FontSize',12);
title('Density of minDt on Membrane','FontSize',16);
c=colorbar;
ylabel(c,'N_{minDt}(t)');

% EMinDt Plot
figure(3)
surf(time,distance,EminDt');
axis([0 1200 0 6]);
set(gca, 'FontSize',12);
xlabel('time(s)', 'FontSize',12);
ylabel('distance along major axis({\mu}m)','FontSize',12);
title('Density of EminDt on Membrane','FontSize',16);
c=colorbar;
ylabel(c,'N_{EminDt}(t)');