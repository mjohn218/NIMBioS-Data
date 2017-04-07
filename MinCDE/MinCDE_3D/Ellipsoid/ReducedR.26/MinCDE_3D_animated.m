%% MinCDE2D Animated Line

%% Import
load('minDt_R26.mat');
load('EminDt_R26.mat');
load('time_R26.mat');
load('distance_R26.mat');

%% Animation
figure(4)
clf
MinDt = animatedline('Color',[0 .75 1],'LineWidth',3);
MinEDt = animatedline('Color',[0 0 1],'LineWidth',3);

axis([0 5.5 0 100])
xlabel('Distance along long axis (um)','FontSize',14);
ylabel('N(x)','FontSize',14);
title('Changes in Number of Particles over Time in MinCDE3D_{R26}','FontSize',16);
l=legend('MinDt','EminDT');
l.FontSize = 14;
for k = 1:length(time_R26)*8
    clearpoints(MinDt)
    clearpoints(MinEDt)
    
    index = uint8(ceil(k/8));
    addpoints(MinDt,distance_R26,minDt_R26(index,:))
    addpoints(MinEDt,distance_R26,EminDt_R26(index,:))
    
    drawnow
end
