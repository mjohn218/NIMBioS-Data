%% Figure 2 - Reversible Reaction 3D Plot
% Clear
close all
clear
clc

%% Import Raw Data
% t is in us
load('time.mat');
load('det1D_1000.mat');
load('det3D_1000.mat');
load('gillespie_1000.mat');
load('mcell_1000.mat');
load('smoldyn_1000.mat');

load('det1D_200.mat');
load('det3D_200.mat');
load('gillespie_200.mat');
load('mcell_200.mat');
load('smoldyn_200.mat');

%% Process data
%   Data files containing one or more trials are averaged
%   Data points corresponding to the time vector are extrapolated
%   Time is converted to us

% ODE_1000
%   Has duplicated values near the end
%   time: s
det1D_1000 = interp1(det1D_1000(:,1)*1e6, det1D_1000(:,2),time,'pchip',0);

% PDE_1000
%   2 trials
%   time: s
det3D_1000(:,2) = (det3D_1000(:,2) + det3D_1000(:,4))/2;
det3D_1000 = interp1(det3D_1000(:,1)*1e6,det3D_1000(:,2), time,'pchip', 0);

% Gillespie_1000
%   20 trials
%   time: s
temp = 0;
for n = 1:2:20
    temp = temp(:,1) + interp1(gillespie_1000(:,n)*1e6,gillespie_1000(:,n+1), time, 'pchip',0);
end
gillespie_1000 = temp(:,1)/10;
clear temp

% MCell_1000
%   time: s

% Smoldyn_1000
smoldyn_1000 = interp1(smoldyn_1000(:,1)*1e6, smoldyn_1000(:,2), time, 'pchip', 0);

% ODE_200
%   time: s
det1D_200 = interp1(det1D_200(:,1)*1e6,det1D_200(:,2), time, 'pchip',0);

% PDE_200
%   2 trials
%   time: s
det3D_200(:,2) = (det3D_200(:,2) + det3D_200(:,4))/2;
det3D_200 = interp1(det3D_200(:,1)*1e6,det3D_200(:,2), time, 'pchip', 0);

% Gillespie_200
%   20 trials
%   time: s
temp = 0;
for n = 1:2:20
    temp = temp(:,1) + interp1(gillespie_200(:,n)*1e6,gillespie_200(:,n+1), time, 'pchip',0);
end
gillespie_200 = temp(:,1)/10;
clear temp

% MCellBC_200
%   time: s


% Smoldyn_200
smoldyn_200 = interp1(smoldyn_200(:,1)*1e6, smoldyn_200(:,2), time, 'pchip', 0);

%% Plot
hold on
% Ka = 1000 nm3.us-1
semilogx(time,det1D_1000, '-','Color',[0 0 0], 'LineWidth', 1);
semilogx(time,det3D_1000, '-','Color',[0 0 1], 'LineWidth', 1);
semilogx(time,gillespie_1000, '-','Color',[0 .75 1], 'LineWidth', 1);
semilogx(time,smoldyn_1000, '-','Color',[0 .8 0], 'LineWidth', 1);
semilogx(time,mcell_1000, '-','Color',[.5 0 .8], 'LineWidth', 1);

% Absorbing BC
semilogx(time,det1D_200, '--','Color',[0 0 0], 'LineWidth', 1);
semilogx(time,det3D_200, '--','Color',[0 0 1], 'LineWidth', 1);
semilogx(time,gillespie_200, '--','Color',[0 .75 1], 'LineWidth', 1);
semilogx(time,smoldyn_200, '--','Color',[0 .8 0], 'LineWidth', 1);
semilogx(time,mcell_200, '--','Color',[.5 0 .8], 'LineWidth', 1);

set(gca, 'xscale', 'log', 'yscale','log');
axis([0 1e8 0 1e5]);
title('Reversible Reaction');
xlabel('Time (us)');
ylabel('N_A(t)');
legend('ODE Det., K_a = 1000 nm^3/\mus', 'PDE Det., K_a = 1000 nm^3/\mus','Gillespie, K_a = 1000 nm^3/\mus', 'Smoldyn, K_a = 1000 nm^3/\mus','MCell, K_a = 1000 nm^3/\mus','ODE Det., K_a = 200 nm^3\mus','PDE Det., K_a = 200 nm^3\mus', 'Gillespie, K_a = 200 nm^3\mus', 'Smoldyn, K_a = 200 nm^3\mus', 'MCell, K_a = 200 nm^3\mus');

% Legend        
%   ODE Deterministic       Black           [0 0 0]
%   PDE Deterministic       Dark Blue       [0 0 1]
%   Gillespie               Light Blue      [0 .75 1]
%   FPR                     Yellow          [1 .85 0]
%   Smoldyn                 Green           [0 .8 0]
%   MCell                   Purple          [.5 0 .8]
