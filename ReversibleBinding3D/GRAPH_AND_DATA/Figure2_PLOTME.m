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
load('fpr_1000.mat');

load('det1D_200.mat');
load('det3D_200.mat');
load('gillespie_200.mat');
load('mcell_200.mat');
load('smoldyn_200.mat');
load('fpr_200.mat');

%% Process data
%   Data files containing one or more trials are averaged
%   Data points corresponding to the time vector are extrapolated
%   Time is converted to us

% ODE_1000
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

% FPR_1000
%   time: us
fpr_1000 = interp1(fpr_1000(:,1),fpr_1000(:,2),time,'pchip',0);

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

% MCell_200
%   time: s

% Smoldyn_200
smoldyn_200 = interp1(smoldyn_200(:,1)*1e6, smoldyn_200(:,2), time, 'pchip', 0);

% FPR_200
fpr_200 = interp1(fpr_200(:,1), fpr_200(:,2),time,'pchip',0);

%% Plot
% Ka = 1000 nm^3.us-1
figure(1)
subplot(3,1,1) % Spatial Effects
hold on
title('Spatial Effects');
set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize',12, 'fontweight','bold',...
    'YTick',[1e1 1e2 1e3 1e4 1e5]);
axis([0 1e6 0 1e5]);
g1 = semilogx(time,det1D_1000, '-','Color',[0 0 1], 'LineWidth', 4);
g2 = semilogx(time,det3D_1000, '--','Color',[0 .75 1], 'LineWidth', 3);

subplot(3,1,2) % Stochastic Effects
hold on
title('Stochastic Effects');
set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize',12, 'fontweight','bold',...
    'YTick',[1e1 1e2 1e3 1e4 1e5]);
axis([0 1e6 0 1e5]);
ylabel('N_A(t)');
g3 = semilogx(time,det1D_1000, '-','Color',[0 0 1], 'LineWidth', 4);
g4 = semilogx(time,gillespie_1000, '--','Color',[0 .8 0], 'LineWidth', 3);

subplot(3,1,3) % Single Particle Methods
hold on
title('Single Particle Methods');
set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize',12, 'fontweight','bold',...
    'YTick',[1e1 1e2 1e3 1e4 1e5]);
axis([0 1e6 0 1e5]);
xlabel('Time (us)');
g5 = semilogx(time,det1D_1000, '-','Color',[0 0 1], 'LineWidth', 4);
g6 = semilogx(time,mcell_1000, '-.','Color',[.5 0 .8], 'LineWidth', 3);
g7 = semilogx(time,smoldyn_1000, '--','Color',[1 .5 0], 'LineWidth', 3);
g8 = semilogx(time,fpr_1000, '--','Color',[1 .85 0], 'LineWidth', 3);

lgnd = legend([g1 g2 g4 g6 g7 g8],'ODE','PDE','Gillespie',...
    'MCell','Smoldyn','FPR');
lgnd.FontSize = 11.5;
lgnd.Position = [.02 .475 1 1];
lgnd.Orientation = 'horizontal';

% Ka = 200 nm^3.us-1
figure(2)
subplot(3,1,1) % Spatial Effects
hold on
title('Spatial Effects');
set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize',12, 'fontweight','bold',...
    'YTick',[1e1 1e2 1e3 1e4 1e5]);
axis([0 1e6 0 1e5]);
g1b = semilogx(time,det1D_200, '-','Color',[0 0 1], 'LineWidth', 4);
g2b = semilogx(time,det3D_200, '--','Color',[0 .75 1], 'LineWidth', 3);

subplot(3,1,2) % Stochastic Effects
hold on
title('Stochastic Effects');
set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize',12, 'fontweight','bold',...
    'YTick',[1e1 1e2 1e3 1e4 1e5]);
axis([0 1e6 0 1e5]);
ylabel('N_A(t)');
g3b = semilogx(time,det1D_200, '-','Color',[0 0 1], 'LineWidth', 4);
g4b = semilogx(time,gillespie_200, '--','Color',[0 .8 0], 'LineWidth', 3);

subplot(3,1,3) % single Particle Methods
hold on
title('Single Particle Methods');
set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize',12, 'fontweight','bold',...
    'YTick',[1e1 1e2 1e3 1e4 1e5]);
axis([0 1e6 0 1e5]);
xlabel('Time (us)');
g5b = semilogx(time,det1D_200, '-','Color',[0 0 1], 'LineWidth', 4);
g6b = semilogx(time,mcell_200, '-.','Color',[.5 0 .8], 'LineWidth', 3);
g7b = semilogx(time,smoldyn_200, '--','Color',[1 .5 0], 'LineWidth', 3);
g8b = semilogx(time,fpr_200, '--','Color',[1 .85 0], 'LineWidth', 3);

lgnd = legend([g1b g2b g4b g6b g7b g8b],'ODE','PDE','Gillespie',...
    'MCell','Smoldyn','FPR');
lgnd.FontSize = 11.5;
lgnd.Position = [.02 .475 1 1];
lgnd.Orientation = 'horizontal';

% Legend 
%   Theory                  Black           [0 0 0]
%   ODE Deterministic       Dark Blue       [0 0 1]
%   PDE Deterministic       Light Blue      [0 .75 1]
%   Gillespie               Green           [0 .8 0]
%   FPR                     Yellow          [1 .85 0]
%   Smoldyn                 Orange          [1 .5 0]
%   MCell                   Purple          [.5 0 .8]
%% Clear
%clear
%clc