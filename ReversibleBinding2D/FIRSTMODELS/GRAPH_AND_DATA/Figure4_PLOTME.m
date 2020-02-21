%% Figure 3 - Simple Irreversible Reaction 2D Plots
% Clear


%% Import
% t is in us
load('time_1s.dat');
%load('theory_1.mat');
load('det1D_1.mat');
load('det3D_1.mat');
load('gillespie_1.mat');
load('mcell_1.mat');
load('smoldyn_1.mat');
load('fpr_1_traj50.dat');

%load('theory_10.mat');
load('det1D_10.mat');
load('det3D_10.mat');
load('gillespie_10.mat');
load('mcell_10.mat');
load('smoldyn_10.mat');
load('fpr_10_traj50.dat');

%% Process data
%   Data files containing one or more trials are averaged
%   Data points corresponding to the time vector are extrapolated
%   Time is converted to us

% Theory_1 USE FPR
%   time:us
theory_1 = fpr_1_traj50;
%interp1(fpr_1_traj50(:,1), fpr_1_traj50(:,2), time,'pchip',0);

time=time_1s;
% ODE_1
%   time: s
%det1D_1 = interp1(det1D_1(:,1)*1e6, det1D_1(:,2),time,'pchip',0);

% PDE_1
% 
%   time: s

%det3D_1 = interp1(det3D_1(:,1)*1e6,det3D_1(:,2), time,'pchip', 0);

% Gillespie_1
%   20 trials
%   time: s
temp = 0;
for n = 1:2:20
    temp = temp(:,1) + interp1(gillespie_1(:,n)*1e6,gillespie_1(:,n+1), time, 'pchip',0);
end
gillespie_1 = temp(:,1)/10;
tot=[time, gillespie_1];
save gillespie_ka1_10trajavg.dat tot -ASCII
clear temp

% MCell_1
%   10 trials, already averaged
%   time: s
%mcell_1 = interp1(mcell_1(:,1)*1e6, mcell_1(:,2), time, 'pchip', 0);

% Smoldyn_1
%   10 trials
%   time: s
for n = 3:2:20
    smoldyn_1(:,2) = smoldyn_1(:,2)+smoldyn_1(:,n+1);
end
smoldyn_1(:,2) = smoldyn_1(:,2)./10;
save smoldyn_ka1_10trajavg.dat smoldyn_1 -ASCII
%smoldyn_1 = interp1(smoldyn_1(:,1)*1e6, smoldyn_1(:,2), time, 'pchip', 0);

% FPR_1
%   10 trials
%   time: s

% Theory_10
%   time: us
theory_10 = fpr_10_traj50;
%theory_10=interp1(fpr_10_traj50(:,1),fpr_10_traj50(:,2),time,'pchip',0);

% ODE_10
%   time: s
%det1D_10 = interp1(det1D_10(:,1)*1e6,det1D_10(:,2), time, 'pchip',0);

% PDE_10
%   2 trials
%   time: s
%det3D_10(:,2) = (det3D_10(:,2) + det3D_10(:,4))/2;
%det3D_10 = interp1(det3D_10(:,1)*1e6,det3D_10(:,2), time, 'pchip', 0);

% Gillespie_10
%   20 trials
%   time: s
temp = 0;
for n = 1:2:20
    temp = temp(:,1) + interp1(gillespie_10(:,n)*1e6,gillespie_10(:,n+1), time, 'pchip',0);
end
gillespie_10 = temp(:,1)/10;
tot=[time, gillespie_10];
save gillespie_ka10_10traj.dat tot -ASCII

clear temp

% MCell_10
%   time: s
%mcell_10 = interp1(mcell_10(:,1)*1e6, mcell_10(:,2), time, 'pchip', 0);

% Smoldyn_10
%   10 trials
%   time: s
for n = 3:2:20
    smoldyn_10(:,2) = smoldyn_10(:,2)+smoldyn_10(:,n+1);
end
smoldyn_10(:,2) = smoldyn_10(:,2)./10;
save smoldyn_ka10_10traj_avg.dat smoldyn_10 -ASCII
%smoldyn_10 = interp1(smoldyn_10(:,1)*1e6, smoldyn_10(:,2), time, 'pchip', 0);



%% Plot
% Ka = 1 um^2.s-1
figure(8)
subplot(3,1,1) % Spatial Effects
% x0=10;
% y0=10;
% width=400;
% height=590;
% set(gcf,'units','points','position',[x0,y0,width,height])
hold on
title('Spatial Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([1e4 1e6 0 50]);
g0 = semilogx(theory_1(:,1),theory_1(:,2),'-','Color',[0 0 0], 'LineWidth',4);
g1 = semilogx(det1D_1(:,1)*1E6,det1D_1(:,2), '-.','Color',[0 0 1], 'LineWidth', 3);
g2 = semilogx(det3D_1(:,1)*1E6,det3D_1(:,2), '--','Color',[0 .75 1], 'LineWidth', 3);

subplot(3,1,2) % Stochastic Effects
hold on
title('Stochastic Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([1e4 1e6 0 50]);
ylabel('N_A(t)');
g3 = semilogx(theory_1(:,1),theory_1(:,2),'-','Color',[0 0 0], 'LineWidth',4);
g4 = semilogx(det1D_1(:,1)*1E6,det1D_1(:,2), '-.','Color',[0 0 1], 'LineWidth', 3);
g5 = semilogx(time,gillespie_1, '--','Color',[0 .8 0], 'LineWidth', 3);

subplot(3,1,3) % Single Particle Methods
hold on
title('Single Particle Methods');
axis([1e4 1e6 0 50]);
xlabel('Time (us)');
g6 = semilogx(theory_1(:,1),theory_1(:,2),'-','Color',[0 0 0], 'LineWidth',4);
g7 = semilogx(mcell_1(:,1)*1E6,mcell_1(:,2), '--','Color',[.5 0 .8], 'LineWidth', 3);
g8 = semilogx(smoldyn_1(:,1)*1E6,smoldyn_1(:,2), '-.','Color',[1 .5 0], 'LineWidth', 3);
%g9 = semilogx(time,fpr_1, '--','Color',[1 .85 0], 'LineWidth', 3);

lgnd = legend([g0 g1 g2 g5 g7 g8],'FPR','ODE','PDE','Gillespie',...
    'MCell', 'Smoldyn');
lgnd.FontSize = 11.5;
lgnd.Position = [.02 .475 1 1];
lgnd.Orientation = 'horizontal';

% Ka = 10 um^2.s-1
figure(7)
subplot(3,1,1) % Spatial Effects
% x0=10;
% y0=10;
% width=400;
% height=590;
% set(gcf,'units','points','position',[x0,y0,width,height])
hold on
title('Spatial Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([1e4 1e6 0 50]);
g0 = semilogx(theory_10(:,1),theory_10(:,2),'-','Color',[0 0 0], 'LineWidth',4);
g1 = semilogx(det1D_10(:,1)*1E6,det1D_10(:,2), '-.','Color',[0 0 1], 'LineWidth', 3);
g2 = semilogx(det3D_10(:,1)*1E6,det3D_10(:,2), '--','Color',[0 .75 1], 'LineWidth', 3);

subplot(3,1,2) % Stochastic Effects
hold on
title('Stochastic Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([1e4 1e6 0 50]);
ylabel('N_A(t)');
g3 = semilogx(theory_10(:,1),theory_10(:,2),'-','Color',[0 0 0], 'LineWidth',4);
g4 = semilogx(det1D_10(:,1)*1E6,det1D_10(:,2), '-','Color',[0 0 1], 'LineWidth', 3);
g5 = semilogx(time,gillespie_10, '--','Color',[0 .8 0], 'LineWidth', 3);

subplot(3,1,3) % Single Particle Methods
hold on
title('Single Particle Methods');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([1e4 1e6 0 50]);
xlabel('Time (us)');
g6 = semilogx(theory_10(:,1),theory_10(:,2),'-','Color',[0 0 0], 'LineWidth',4);
g7 = semilogx(mcell_10(:,1)*1E6,mcell_10(:,2), '--','Color',[.5 0 .8], 'LineWidth', 3);
g8 = semilogx(smoldyn_10(:,1)*1E6,smoldyn_10(:,2), '-.','Color',[1 .5 0], 'LineWidth', 3);
%g9 = semilogx(time,fpr_10, '--','Color',[1 .85 0], 'LineWidth', 3);

lgnd = legend([g0 g1 g2 g5 g7 g8 ],'FPR','ODE','PDE','Gillespie',...
    'MCell','Smoldyn');
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