%% Figure 3 - Simple Irreversible Reaction 2D Plots
% Clear


%% Import
load('theory_1.mat'); % t in us
load('det1D_1.mat'); % t in s
load('det3D_1.mat'); % t in s
load('gillespie_1.mat'); % t in s
load('mcell_1.mat'); % t in s
load('smoldyn_1.mat'); % t in s
load('fpr_1.mat'); % t in s

load('theory_10.mat'); % t in us
load('det1D_10.mat'); % t in s
load('det3D_10.mat'); % t in s
load('gillespie_10.mat'); % t in s
load('mcell_10.mat'); % t in s
load('smoldyn_10.mat'); % t in s
load('fpr_10.mat'); % t in s

load('theory_100.mat'); % t in us
load('det1D_100.mat'); % t in s
load('det3D_100.mat'); % t in s
load('gillespie_100.mat'); % t in s
load('mcell_100.mat'); % t in s
load('smoldyn_100.mat'); % t in s
load('fpr_100.mat'); % t in s
%% Process data
%   Data files containing one or more trials are averaged
%   Time is converted to us

% ODE_1
%   time: s
det1D_1(:,1) = det1D_1(:,1)*1e6;

% PDE_1
%   2 trials
%   time: s
det3D_1(:,2) = (det3D_1(:,2) + det3D_1(:,4))/2;
det3D_1(:,1) = det3D_1(:,1)*1e6;

% Gillespie_1
%   20 trials
%   time: s
temp = 0;
for n = 1:2:20
    temp = temp(:,1) + interp1(gillespie_1(:,n)*1e6,gillespie_1(:,n+1), det1D_1(:,1), 'pchip',0);
end
clear gillespie_1
gillespie_1(:,1) = det1D_1(:,1);
gillespie_1(:,2) = temp(:,1)/10;
clear temp

% MCell_1
%   time: s
mcell_1(:,1) = mcell_1(:,1)*1e6;

% Smoldyn_1
%   10 trials
%   time: s
for n = 3:2:20
    smoldyn_1(:,2) = smoldyn_1(:,2)+smoldyn_1(:,n+1);
end
smoldyn_1(:,2) = smoldyn_1(:,2)./10;
smoldyn_1(:,1)= smoldyn_1(:,1)*1e6;

% FPR_1
%   24 trials
%   time: s
for n = 3:25
    fpr_1(:,2) = fpr_1(:,2) + fpr_1(:,n);
end
fpr_1(:,2) = fpr_1(:,2)./24;
fpr_1(:,1) = fpr_1(:,1)*1e6;

% ODE_10
%   time: s
det1D_10(:,1) = det1D_10(:,1)*1e6;

% PDE_10
%   2 trials
%   time: s
det3D_10(:,2) = (det3D_10(:,2) + det3D_10(:,4))/2;
det3D_10(:,1) = det3D_10(:,1)*1e6;

% Gillespie_10
%   20 trials
%   time: s
temp = 0;
for n = 1:2:20
    temp = temp(:,1) + interp1(gillespie_10(:,n)*1e6,gillespie_10(:,n+1), det1D_10(:,1), 'pchip',0);
end
clear gillespie_10
gillespie_10(:,2) = det1D_10(:,1);
gillespie_10(:,2) = temp(:,1)/10;
clear temp

% MCell_10
%   time: s
mcell_10(:,1) = mcell_10(:,1)*1e6;

% Smoldyn_10
%   10 trials
%   time: s
for n = 3:2:20
    smoldyn_10(:,2) = smoldyn_10(:,2)+smoldyn_10(:,n+1);
end
smoldyn_10(:,2) = smoldyn_10(:,2)./10;
smoldyn_10(:,1) = smoldyn_10(:,1)*1e6;

% FPR_10
%   24 trials
%   time: s
for n = 3:25
    fpr_10(:,2) = fpr_10(:,2) + fpr_10(:,n);
end
fpr_10(:,2) = fpr_10(:,2)./24;
fpr_10(:,1) = fpr_10(:,1)*1e6;

% ODE_100
%   time: s
det1D_100(:,1) = det1D_100(:,1)*1e6;

% PDE_100
%   2 trials
%   time: s
det3D_100(:,2) = (det3D_100(:,2) + det3D_100(:,4))/2;
det3D_100(:,1) = det3D_100(:,1)*1e6;

% Gillespie_100
%   20 trials
%   time: s
temp = 0;
for n = 1:2:20
    temp = temp(:,1) + interp1(gillespie_100(:,n)*1e6,gillespie_100(:,n+1), det1D_100(:,1), 'pchip',0);
end
clear gillespie_100;
gillespie_100(:,1) = det1D_100(:,1);
gillespie_100(:,2) = temp(:,1)/10;
clear temp

% MCell_100
%   time: s
mcell_100(:,1)= mcell_100(:,1)*1e6;

% Smoldyn_100
%   10 trials
%   time: s
for n = 3:2:20
    smoldyn_100(:,2) = smoldyn_100(:,2)+smoldyn_100(:,n+1);
end
smoldyn_100(:,2) = smoldyn_100(:,2)./10;
smoldyn_100(:,1) = smoldyn_100(:,1)*1e6;

% FPR_100
%   Deletes first row
%   10 trials
%   time: s
for n = 3:25
    fpr_100(:,2) = fpr_100(:,2) + fpr_100(:,n);
end
fpr_100(:,2) = fpr_100(:,2)./24;
fpr_100(:,1) = fpr_100(:,1)*1e6;
%% Plot
% Ka = 1 um^2.s-1
figure(3)
subplot(3,1,1) % Spatial Effects
hold on
title('Spatial Effects');
set(gca,'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 1.5e5 0 1000]);
g0 = semilogx(theory_1(:,1),theory_1(:,2),'-','Color',[0 0 0], 'LineWidth', 4);
g1 = semilogx(det1D_1(:,1),det1D_1(:,2), '-.','Color',[0 0 1], 'LineWidth', 3);
g2 = semilogx(det3D_1(:,1),det3D_1(:,2), '--','Color',[0 .75 1], 'LineWidth', 3);

subplot(3,1,2) % Stochastic Effects
hold on
title('Stochastic Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 1e5 0 1000]);
ylabel('N_A(t)');
g3 = semilogx(theory_1(:,1),theory_1(:,2),'-','Color',[0 0 0], 'LineWidth', 4);
g4 = semilogx(det1D_1(:,1),det1D_1(:,2), '-.','Color',[0 0 1], 'LineWidth', 3);
g5 = semilogx(gillespie_1(:,1),gillespie_1(:,2), '--','Color',[0 .8 0], 'LineWidth', 3);

subplot(3,1,3) % Single Particle Methods
hold on
title('Single Particle Methods');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 1e5 0 1000]);
xlabel('Time (us)');
g6 = semilogx(theory_1(:,1),theory_1(:,2),'-','Color',[0 0 0], 'LineWidth', 4);
g7 = semilogx(mcell_1(:,1),mcell_1(:,2), '--','Color',[.5 0 .8], 'LineWidth', 3);
g8 = semilogx(smoldyn_1(:,1),smoldyn_1(:,2), '--','Color',[1 .5 0], 'LineWidth', 3);
g9 = semilogx(fpr_1(:,1),fpr_1(:,2), '--','Color',[1 .85 0], 'LineWidth', 3);

lgnd = legend([g0 g1 g2 g5 g7 g8 g9],'Theory','ODE','PDE','Gillespie',...
    'MCell', 'Smoldyn', 'FPR');
lgnd.FontSize = 11.5;
lgnd.Position = [.02 .475 1 1];
lgnd.Orientation = 'horizontal';

% Ka = 10 um^2.s-1
figure(4)
subplot(3,1,1) % Spatial Effects
% x0=10;
% y0=10;
% width=400;
% height=590;
% set(gcf,'units','points','position',[x0,y0,width,height])
hold on
title('Spatial Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 1e5 0 1000]);
g0 = semilogx(theory_10(:,1),theory_10(:,2),'-','Color',[0 0 0], 'LineWidth', 4);
g1 = semilogx(det1D_10(:,1),det1D_10(:,2), '-.','Color',[0 0 1], 'LineWidth', 3);
g2 = semilogx(det3D_10(:,1),det3D_10(:,2), '--','Color',[0 .75 1], 'LineWidth', 3);

subplot(3,1,2) % Stochastic Effects
hold on
title('Stochastic Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 1e5 0 1000]);
ylabel('N_A(t)');
g3 = semilogx(theory_10(:,1),theory_10(:,2),'-','Color',[0 0 0], 'LineWidth', 4);
g4 = semilogx(det1D_10(:,1),det1D_10(:,2), '-.','Color',[0 0 1], 'LineWidth', 3);
g5 = semilogx(gillespie_10(:,1),gillespie_10(:,2), '--','Color',[0 .8 0], 'LineWidth', 3);

subplot(3,1,3) % Single Particle Methods
hold on
title('Single Particle Methods');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 1e5 0 1000]);
xlabel('Time (us)');
g6 = semilogx(theory_10(:,1),theory_10(:,2),'-','Color',[0 0 0], 'LineWidth', 4);
g7 = semilogx(mcell_10(:,1),mcell_10(:,2), '--','Color',[.5 0 .8], 'LineWidth', 3);
g8 = semilogx(smoldyn_10(:,1),smoldyn_10(:,2), '-.','Color',[1 .5 0], 'LineWidth', 3);
g9 = semilogx(fpr_10(:,1),fpr_10(:,2), '--','Color',[1 .85 0], 'LineWidth', 3);

lgnd = legend([g0 g1 g2 g5 g7 g8 g9],'Theory','ODE','PDE','Gillespie',...
    'MCell', 'Smoldyn','FPR');
lgnd.FontSize = 11.5;
lgnd.Position = [.02 .475 1 1];
lgnd.Orientation = 'horizontal';

% Ka = 100 um^2.s-1
figure(5)
subplot(3,1,1) % Spatial Effects
% x0=10;
% y0=10;
% width=400;
% height=590;
% set(gcf,'units','points','position',[x0,y0,width,height])
hold on
title('Spatial Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 1e6 0 1000]);
g0 = semilogx(theory_100(:,1),theory_100(:,2),'-','Color',[0 0 0], 'LineWidth', 4);
g1 = semilogx(det1D_100(:,1),det1D_100(:,2), '-.','Color',[0 0 1], 'LineWidth', 3);
g2 = semilogx(det3D_100(:,1),det3D_100(:,2), '--','Color',[0 .75 1], 'LineWidth', 3);

subplot(3,1,2) % Stochastic Effects
hold on
title('Stochastic Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 1e5 0 1000]);
ylabel('N_A(t)');
g3 = semilogx(theory_100(:,1),theory_100(:,2),'-','Color',[0 0 0], 'LineWidth', 4);
g4 = semilogx(det1D_100(:,1),det1D_100(:,2), '-.','Color',[0 0 1], 'LineWidth', 3);
g5 = semilogx(gillespie_100(:,1),gillespie_100(:,2), '--','Color',[0 .8 0], 'LineWidth', 3);

subplot(3,1,3) % Single Particle Methods
hold on
title('Single Particle Methods');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([0 1e5 0 1000]);
xlabel('Time (us)');
g6 = semilogx(theory_100(:,1),theory_100(:,2),'-','Color',[0 0 0], 'LineWidth', 4);
g7 = semilogx(mcell_100(:,1),mcell_100(:,2), '--','Color',[.5 0 .8], 'LineWidth', 3);
g8 = semilogx(smoldyn_100(:,1),smoldyn_100(:,2), '-.','Color',[1 .5 0], 'LineWidth', 3);
g9 = semilogx(fpr_100(:,1),fpr_100(:,2), '--','Color',[1 .85 0], 'LineWidth', 3);

lgnd = legend([g0 g1 g2 g5 g7 g8 g9],'Theory','ODE','PDE','Gillespie',...
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