%% PLOTME_3DIrrBinding_DiffConc.m
%  Plotting script for the 3D irreversible binding - diffusion controlled
%  parameter set. 
%  Use: PLOTME_3DIrrBinding_DiffConc.m

%% Initialize



%% Import
load ODE
load Gillespie
load PDE
load Smoldyn
load MCell
load FPR

%% Data Processing
%  Average stochastic simulations
%  Convert time to us from s

% Theory
%   time: us
time = MCell(:,1)*1e6;
initA = 1e-2; %initial concentration in nm^-3
sigma = 1; %nm
D = 100; %nm^2.um-1
kD = 4*pi*sigma*D;
ka = 1000; %nm^3.um-1
lambda = (D/sigma^2)*(1+(ka/kD))^2;
Theory = ((1/initA)-kD*(ka/(ka+kD))*(ka/(kD*lambda))+kD*(ka/(ka+kD)).*...
    (time+((ka.*2.*sqrt(time))./(kD.*sqrt(lambda*pi)))+...
    (ka/(kD*lambda)).*erfcx(sqrt(lambda*time)))).^(-1);
Theory = [time Theory*2e5]; % Convert to number of molecules
clear('time','initA','sigma','D','kD','ka','lambda');

% ODE
ODE(:,1) = ODE(:,1)*1e6;

% Gillespie
tempt = zeros(998,1);
for r = 1:10
    tempt = tempt + Gillespie(1:998,2*r - 1);
end
tempt = tempt/10;
Gillespie = [tempt*1e6 Gillespie(1:998,2)];
clear tempt

% PDE
PDE(:,1) = PDE(:,1)*1e6;

% Smoldyn
temp = zeros(1100,1);
for r = 1:10
    temp = temp + Smoldyn(:,2*r);
end
temp = temp/10;
Smoldyn = [Smoldyn(:,1)*1e6 temp];
clear temp

% MCell
temp = zeros(length(MCell),1);
for r = 1:10
    temp = temp + MCell(:,2*r);
end
temp = temp/10;
MCell = [MCell(:,1)*1e6 temp];
clear temp

% FPR
%   Deletes first row
%   10 trials, already averaged
%   time: us
FPR = FPR(2:length(FPR),:);

%% Error Calculation
% Theory Initialization
initA = 1e-2; %initial concentration in nm^-3
sigma = 1; %nm
D = 100; %nm^2.um-1
kD = 4*pi*sigma*D;
ka = 1000; %nm^3.um-1
lambda = (D/sigma^2)*(1+(ka/kD))^2;

% ODE
time = ODE(:,1);
exact = ((1/initA)-kD*(ka/(ka+kD))*(ka/(kD*lambda))+kD*(ka/(ka+kD)).*...
    (time+((ka.*2.*sqrt(time))./(kD.*sqrt(lambda*pi)))+...
    (ka/(kD*lambda)).*erfcx(sqrt(lambda*time)))).^(-1);
exact = exact*2e5; % Convert to number of molecules
ODEerr = abs(ODE(:,2)-exact)./exact;

% Gillespie
time = Gillespie(:,1);
exact = ((1/initA)-kD*(ka/(ka+kD))*(ka/(kD*lambda))+kD*(ka/(ka+kD)).*...
    (time+((ka.*2.*sqrt(time))./(kD.*sqrt(lambda*pi)))+...
    (ka/(kD*lambda)).*erfcx(sqrt(lambda*time)))).^(-1);
exact = exact*2e5; % Convert to number of molecules
Gillespieerr = abs(Gillespie(:,2)-exact)./exact;

% PDE
time = PDE(:,1);
exact = ((1/initA)-kD*(ka/(ka+kD))*(ka/(kD*lambda))+kD*(ka/(ka+kD)).*...
    (time+((ka.*2.*sqrt(time))./(kD.*sqrt(lambda*pi)))+...
    (ka/(kD*lambda)).*erfcx(sqrt(lambda*time)))).^(-1);
exact = exact*2e5; % Convert to number of molecules
PDEerr = abs(PDE(:,2)-exact)./exact;

% Smoldyn
time = Smoldyn(:,1);
exact = ((1/initA)-kD*(ka/(ka+kD))*(ka/(kD*lambda))+kD*(ka/(ka+kD)).*...
    (time+((ka.*2.*sqrt(time))./(kD.*sqrt(lambda*pi)))+...
    (ka/(kD*lambda)).*erfcx(sqrt(lambda*time)))).^(-1);
exact = exact*2e5; % Convert to number of molecules
Smoldynerr = abs(Smoldyn(:,2)-exact)./exact;

% MCell
time = MCell(:,1);
exact = ((1/initA)-kD*(ka/(ka+kD))*(ka/(kD*lambda))+kD*(ka/(ka+kD)).*...
    (time+((ka.*2.*sqrt(time))./(kD.*sqrt(lambda*pi)))+...
    (ka/(kD*lambda)).*erfcx(sqrt(lambda*time)))).^(-1);
exact = exact*2e5; % Convert to number of molecules
MCellerr = abs(MCell(:,2)-exact)./exact;

% FPR
time = FPR(:,1);
exact = ((1/initA)-kD*(ka/(ka+kD))*(ka/(kD*lambda))+kD*(ka/(ka+kD)).*...
    (time+((ka.*2.*sqrt(time))./(kD.*sqrt(lambda*pi)))+...
    (ka/(kD*lambda)).*erfcx(sqrt(lambda*time)))).^(-1);

tot1=[time, exact];
save theory_At.dat tot1 -ASCII
exact = exact*2e5; % Convert to number of molecules
FPRerr = abs(FPR(:,2)-exact)./exact;

tot=[time, exact];
save theory_irr.dat tot -ASCII

% PRINT
algorithm = {'ODE';'PDE';'Gillespie';'Smoldyn';'MCell';'FPR'};
re = {'t1eneg1','t1e0','t1e1','t1e2','eqerr'};
% Average Relative Error from 0s to 1e-1 us
re1 = [rsum(ODEerr,ODE(:,1),0,1e-1);...
    rsum(PDEerr,PDE(:,1),0,1e-1);...
    rsum(Gillespieerr,Gillespie(:,1),0,1e-1);...
    rsum(Smoldynerr,Smoldyn(:,1),0,1e-1);...
    rsum(MCellerr,MCell(:,1),0,1e-1);...
    rsum(FPRerr,FPR(:,1),0,1e-1)];
% Average Relative Error from 1e-1 us to 1e0 us
re2 = [rsum(ODEerr,ODE(:,1),1e-1,1e0);...
    rsum(PDEerr,PDE(:,1),1e-1,1e0);...
    rsum(Gillespieerr,Gillespie(:,1),1e-1,1e0);...
    rsum(Smoldynerr,Smoldyn(:,1),1e-1,1e0);...
    rsum(MCellerr,MCell(:,1),1e-1,1e0);...
    rsum(FPRerr,FPR(:,1),1e-1,1e0)];
% Average Relative Error from 1e0 us to 1e1 us
re3 = [rsum(ODEerr,ODE(:,1),1e0,1e1);
    rsum(PDEerr,PDE(:,1),1e0,1e1);...
    rsum(Gillespieerr,Gillespie(:,1),1e0,1e1);...
    rsum(Smoldynerr,Smoldyn(:,1),1e0,1e1);...
    rsum(MCellerr,MCell(:,1),1e0,1e1);...
    rsum(FPRerr,FPR(:,1),1e0,1e1)];
% Average Relative Error from 1e1 to 1e2 us
re4 = [rsum(ODEerr,ODE(:,1),1e1,1e2);
    rsum(PDEerr,PDE(:,1),1e1,1e2);...
    rsum(Gillespieerr,Gillespie(:,1),1e1,1e2);...
    rsum(Smoldynerr,Smoldyn(:,1),1e1,1e2);...
    rsum(MCellerr,MCell(:,1),1e1,1e2);...
    rsum(FPRerr,FPR(:,1),1e1,1e2)];
% Average Equilibrium Error from 1e2 to end 
eq = [rsum(ODEerr,ODE(:,1),1e2);
    rsum(PDEerr,PDE(:,1),1e2);...
    rsum(Gillespieerr,Gillespie(:,1),1e2);...
    rsum(Smoldynerr,Smoldyn(:,1),1e2);...
    rsum(MCellerr,MCell(:,1),1e2);...
    rsum(FPRerr,FPR(:,1),1e2)];

T = table(re1,re2,re3,re4,eq,'RowNames',algorithm,'VariableNames',re)

clear('time','initA','sigma','D','kD','ka','lambda','re1','re2','re3','re4','re','algorithm');


%% Plot
% Legend 
%   Theory                  Black           [0 0 0]
%   ODE Deterministic       Dark Blue       [0 0 1]
%   PDE Deterministic       Light Blue      [0 .75 1]
%   Gillespie               Green           [0 .8 0]
%   FPR                     Yellow          [1 .85 0]
%   Smoldyn                 Orange          [1 .5 0]
%   MCell                   Purple          [.5 0 .8]

figure(1)
subplot(3,1,1) % Spatial Effects
hold on
title('Spatial Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([1e-3 1e2 0 2000]);
g1 = semilogx(Theory(:,1),Theory(:,2), '-','Color',[0 0 0], 'LineWidth', 4);
g2 = semilogx(ODE(:,1),ODE(:,2), '-.','Color',[0 0 1], 'LineWidth', 4);
g3 = semilogx(PDE(:,1),PDE(:,2), '--','Color',[0 .75 1], 'LineWidth', 4);

subplot(3,1,2) % Stochastic Effects
hold on
title('Stochastic Effects');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([1e-3 1e2 0 2000]);
ylabel('N_A(t)');
g4 = semilogx(Theory(:,1),Theory(:,2), '-','Color',[0 0 0], 'LineWidth', 4);
g5 = semilogx(Gillespie(:,1),Gillespie(:,2), '-.','Color',[0 .8 0], 'LineWidth', 3);
g6 = semilogx(ODE(:,1),ODE(:,2), '--','Color',[0 0 1], 'LineWidth', 3);

subplot(3,1,3) % Single Particle Methods
hold on
title('Single Particle Methods');
set(gca, 'xscale', 'log', 'fontsize',12, 'fontweight','bold');
axis([1e-3 1e2 0 2000]);
xlabel('Time (us)');
g7 = semilogx(Theory(:,1),Theory(:,2), '-','Color',[0 0 0], 'LineWidth', 4);
g8 = semilogx(MCell(:,1),MCell(:,2), '--','Color',[.5 0 .8], 'LineWidth', 3);
g9 = semilogx(Smoldyn(:,1),Smoldyn(:,2), '-.','Color',[1 .5 0], 'LineWidth', 3);
g10 = semilogx(FPR(:,1),FPR(:,2), '--','Color',[1 .85 0], 'LineWidth', 3);

lgnd = legend([g1 g6 g3 g5 g8 g9 g10],'Theory','ODE','PDE','Gillespie',...
    'MCell','Smoldyn','FPR');
lgnd.FontSize = 11.5;
lgnd.Position = [.02 .475 1 1];
lgnd.Orientation = 'horizontal';

%% Reimann Sum Function
% Assumes data is in a column vector
function [reisum] = rsum(data,time,startt,endt)
    % Check number of inputs.
    % Fill in unset optional values.
    switch nargin
        case 2
            startt = time(1);
            endt = time(length(time));
        case 3
            endt = time(length(time));
    end
    
    % Find time window; assign to vector called dist
    dist = time(find(time<=startt,1,'last'):find(time<=endt,1,'last'));
    
    % Calculate the difference in time and the total time covered
    diff = dist - [0; dist(1:size(dist,1)-1)];
    total = dist(size(dist,1)) - dist(1);
    
    % Preallocate Reimann Sum Component Vector
    reimann = zeros(size(dist,1),1);
    % Initialize last value
    reimann(length(reimann)) = data(length(reimann))*1/2*diff(length(reimann));
    
    for i = 1:length(reimann)-1
        reimann(i) = (1/2*diff(i)+1/2*diff(i+1))*data(i);
    end
    
    reisum = sum(reimann)/total;
end

