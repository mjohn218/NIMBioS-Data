%% MinCDE2D Animated Line

%% Import
load('minDt_redL.mat');
load('EminDT_redL.mat');
load('time_redL.mat');
load('distance_redL.mat');

%% Animation
figure(4)
clf
MinDt = animatedline('Color',[0 .75 1],'LineWidth',3);
MinEDt = animatedline('Color',[0 0 1],'LineWidth',3);

axis([0 3 0 200])
xlabel('Distance along long axis (um)','FontSize',14);
ylabel('N(x)','FontSize',14);
title('Changes in Number of Particles over Time in MinCDE3D_{redL}','FontSize',16);
l=legend('MinDt','EminDT');
l.FontSize = 14;
for k = 1:length(time_redL)*8
    clearpoints(MinDt)
    clearpoints(MinEDt)
    
    index = uint8(ceil(k/8));
    addpoints(MinDt,distance_redL,minDt_redL(index,:))
    addpoints(MinEDt,distance_redL,EminDT_redL(index,:))
    
    drawnow
end
