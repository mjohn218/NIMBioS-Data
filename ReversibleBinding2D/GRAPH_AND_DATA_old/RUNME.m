close all
clear
clc
load smoldyn_1
uiopen('C:\Users\osman\SkyDrive (2)\Documents\Projects\CURRENT\paper2\VCell\rdfig_nimbios\JCPFig8\Yogurtcu_2DRD_Fig8_Reversible_withKonRho_Gillespie.fig',1)
hold on

numtraj = 10;
t = [0 10.^(linspace(-6,-1,100))];
At = 0;

for i = 1:numtraj

    At = At + interp1(smoldyn(:,2*(i-1)+1),smoldyn(:,2*(i-1)+2),t);
 
end

At = At/numtraj;
semilogx(t*1e6,At,'k','LineWidth',6);

% 
% % % % % % 
% load det3D_1
% semilogx(t*1e6,At);
% numtraj = 2;
% 
% At = 0;
% 
% for i = 1:numtraj
% 
%     At = At + interp1(det3D(:,2*(i-1)+1),det3D(:,2*(i-1)+2),t);
%  
% end
% 
% At = At/numtraj;
% semilogx(t*1e6,At,'LineWidth',6);
% 
% % % % % % 
% load det1D_1
% semilogx(det1D(:,1)*1e6,det1D(:,2),'LineWidth',6);
% 
% % % % % % 
% load smoldyn_1
% semilogx(smoldyn(:,1)*1e6,smoldyn(:,2),'LineWidth',6);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % CURVE 2
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% load det1D_2
% semilogx(det1D(:,1)*1e6,det1D(:,2),'LineWidth',6);
% 
% % % % % % 
% load det3D_2
% semilogx(t*1e6,At);
% numtraj = 2;
% 
% At = 0;
% 
% for i = 1:numtraj
% 
%     At = At + interp1(det3D(:,2*(i-1)+1),det3D(:,2*(i-1)+2),t);
%  
% end
% 
% At = At/numtraj;
% semilogx(t*1e6,At,'LineWidth',6);
% 
% % % % % % 
load smoldyn_2

numtraj = 10;
t = [0 10.^(linspace(-6,-1,100))];
At = 0;

for i = 1:numtraj

    At = At + interp1(smoldyn(:,2*(i-1)+1),smoldyn(:,2*(i-1)+2),t);
 
end

At = At/numtraj;
semilogx(t*1e6,At,'b','LineWidth',6);

% % % % % % 
% load smoldyn_2
% semilogx(smoldyn(:,1)*1e6,smoldyn(:,2),'LineWidth',6);

saveas(gcf, 'VCELL_MEJ_PRXFig11.fig', 'fig')


