%% Figure 9a - Irreversible Membrane Recruitment
% Clear
close all
clear
clc

%% Import
% t is in us
load('time.mat');
%load('theory_1.mat');
load('det1D.mat');
load('det3D.mat');
load('gillespie.mat');
load('mcell.mat');
load('smoldyn.mat');
load('fpr.mat');

%% Process data
%   Data files containing one or more trials are averaged
%   Data points corresponding to the time vector are extrapolated
%   Time is converted to us

% Theory
%   time:us
%theory = theory(:,2);

% ODE
%   time: s
det1D = interp1(det1D(:,1)*1e6, det1D(:,2),time,'pchip',0);

% PDE
%   2 trials
%   time: s
det3D(:,2) = (det3D(:,2) + det3D(:,4))/2;
det3D = interp1(det3D(:,1)*1e6,det3D(:,2), time,'pchip', 0);

% Gillespie
%   20 trials
%   time: s
temp = 0;
for n = 1:2:20
    temp = temp(:,1) + interp1(gillespie(:,n)*1e6,gillespie(:,n+1), time, 'pchip',0);
end
gillespie = temp(:,1)/10;
clear temp

% MCell
%   1 trial
%   time: s
mcell = interp1(mcell(:,1)*1e6, mcell(:,2), time, 'pchip', 0);

% Smoldyn
%   9 trials
%   time: s
for n = 3:2:18
    smoldyn(:,2) = smoldyn(:,2)+smoldyn(:,n+1);
end
smoldyn(:,2) = smoldyn(:,2)./9;
smoldyn = interp1(smoldyn(:,1)*1e6, smoldyn(:,2), time, 'pchip', 0);

% FPR
%   24 trials, already averaged
%   time: s
fpr = interp1(fpr(:,1)*1e6,fpr(:,2),time,'pchip',0);

%% Plot
% Ka = 1000 um^2.s-1
figure(1)
subplot(3,1,1) % Spatial Effects
% x0=10;
% y0=10;
% width=400;
% height=590;
% set(gcf,'units','points','position',[x0,y0,width,height])
hold on
title('Spatial Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 10000 0 602]);
%g0 = semilogx(time,theory,'-','Color',[0 0 0], 'LineWidth',4);
g2 = semilogx(time,det3D, '-','Color',[0 .75 1], 'LineWidth', 3);
g1 = semilogx(time,det1D, '--','Color',[0 0 1], 'LineWidth', 3);

subplot(3,1,2) % Stochastic Effects
hold on
title('Stochastic Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 10000 0 602]);
ylabel('N_A(t)');
g3 = semilogx(time,det3D,'-','Color',[0 .75 1], 'LineWidth',4);
g4 = semilogx(time,det1D, '-.','Color',[0 0 1], 'LineWidth', 3);
g5 = semilogx(time,gillespie, '--','Color',[0 .8 0], 'LineWidth', 3);

subplot(3,1,3) % Single Particle Methods
hold on
title('Single Particle Methods');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 10000 0 602]);
xlabel('Time (us)');
g6 = semilogx(time,det3D,'-','Color',[0 .75 1], 'LineWidth',4);
g7 = semilogx(time,mcell, '-.','Color',[.5 0 .8], 'LineWidth', 3);
g8 = semilogx(time,smoldyn, ':','Color',[1 .5 0], 'LineWidth', 3);
g9 = semilogx(time,fpr, '--','Color',[1 .85 0], 'LineWidth', 3);

lgnd = legend([g1 g2 g5 g7 g8 g9],'ODE','PDE','Gillespie',...
    'MCell', 'Smoldyn', 'FPR');
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