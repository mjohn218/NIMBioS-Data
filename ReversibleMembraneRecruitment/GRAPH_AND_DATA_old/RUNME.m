close all
clear
clc
load gibsonres_1
% uiopen('C:\Users\osman\SkyDrive (2)\Documents\Projects\CURRENT\paper2\VCell\rdfig_nimbios\PRXFig11\NOINSETPRX_2014_Fig11_Reversible_RDandGillespieWork.fig',1)
% hold on

numtraj = 10;
t = [0 10.^(linspace(-6,log10(4.5),400))];
At = 0;

for i = 1:numtraj

    At = At + interp1(gibsonres(:,2*(i-1)+1),gibsonres(:,2*(i-1)+2),t);
 
end

At = At/numtraj;
semilogx(t*1e6,At,'LineWidth',3);
hold on

% % % % % 
load det3D_cube
% semilogx(t*1e6,At);
numtraj = 2;

t = [0 10.^(linspace(-6,log10(4.5),400))];
At = 0;

for i = 1:numtraj

    At = At + interp1(det3D(:,2*(i-1)+1),det3D(:,2*(i-1)+2),t);
 
end

At = At/numtraj;
semilogx(t*1e6,At,'k','LineWidth',3);

% % % % % % 
% load det3D_1kfover2
% % semilogx(t*1e6,At);
% numtraj = 2;
% 
% t = [0 10.^(linspace(-6,log10(5),100))];
% At = 0;
% 
% for i = 1:numtraj
% 
%     At = At + interp1(det3D(:,2*(i-1)+1),det3D(:,2*(i-1)+2),t);
%  
% end
% 
% At = At/numtraj;
% semilogx(t*1e6,At,'k--','LineWidth',3);

% % % % % 
load det1D_1
semilogx(det1D(:,1)*1e6,det1D(:,2),'LineWidth',3);
% semilogx(t*1e6,interp1(det1D(:,1),det1D(:,2),t),'LineWidth',3);

% % % % % % 
load smoldyn_1
numtraj = 4;

t = [0 10.^(linspace(-6,log10(5),400))];
At = 0;

for i = 1:numtraj

    At = At + interp1(smoldyn(:,2*(i-1)+1),smoldyn(:,2*(i-1)+2),t);
 
end

At = At/numtraj;
semilogx(t*1e6,At,'LineWidth',3);


% % % % % % % % % % % % 

xlabel('Time ({\mu}s)');
ylabel('A(t)')

load A.dat
semilogx(A(:,1)*1e6,A(:,2),'LineWidth',3)
% load FPRoldold.dat
% semilogx(FPRoldold(:,1)*1e6,FPRoldold(:,2)/581*602,'LineWidth',3)
% 
xlim([1 5e6])
legend('ODE Stochastic', 'PDE Deterministic','ODE Deterministic','Smoldyn', 'FPR k_a:1nm^2/{\mu}s','Location','SouthWest')
% title('1 {\mu}M A binds to B (6045 {\mu}m^{-2}) on membrane in 1{\mu}m^3 volume sphere')
ylim([0 1001])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % CURVE 2
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

load det1D_2
semilogx(det1D(:,1)*1e6,det1D(:,2),'LineWidth',6);

% % % % % 
load det3D_2
semilogx(t*1e6,At);
numtraj = 2;

At = 0;

for i = 1:numtraj

    At = At + interp1(det3D(:,2*(i-1)+1),det3D(:,2*(i-1)+2),t);
 
end

At = At/numtraj;
semilogx(t*1e6,At,'LineWidth',6);

% % % % % 
load gibsonres_2

numtraj = 10;
t = [0 10.^(linspace(-6,0,100))];
At = 0;

for i = 1:numtraj

    At = At + interp1(gibsonres(:,2*(i-1)+1),gibsonres(:,2*(i-1)+2),t);
 
end

At = At/numtraj;
semilogx(t*1e6,At,'LineWidth',6);

% % % % % 
load smoldyn_2
semilogx(smoldyn(:,1)*1e6,smoldyn(:,2),'LineWidth',6);

saveas(gcf, 'VCELL_MEJ_PRXFig11.fig', 'fig')


