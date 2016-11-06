clear
clc

load('EminDT.mat');
load('minDt.mat');
load('time.mat');
load('distance.mat');

figure(1)
surf(time,distance,EminDT');
xlabel('time');
ylabel('distance along major axis');
title('Density of EminDT on Membrane');
colorbar

figure(2)
surf(time,distance,minDt');
xlabel('time');
ylabel('distance along major axis');
title('Density of minDt on Membrane');
colorbar