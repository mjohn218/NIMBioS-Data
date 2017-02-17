%% MinCDE2D Animated Line

%% Import
load('minDt.mat');
load('EminDT.mat');
load('time.mat');
load('distance.mat');

%% Animation
figure(4)
clf
MinDt = animatedline('Color',[0 .75 1],'LineWidth',3);
MinEDt = animatedline('Color',[0 0 1],'LineWidth',3);

axis([0 6 0 450])
xlabel('Distance along long axis (um)','FontSize',14);
ylabel('N(x)','FontSize',14);
title('Changes in Number of Particles over Time in MinCDE2D','FontSize',16);
l=legend('MinDt','EminDT');
l.FontSize = 14;
for k = 1:length(time)*8
    clearpoints(MinDt)
    clearpoints(MinEDt)
    
    index = uint8(ceil(k/8));
    addpoints(MinDt,distance,minDt(index,:))
    addpoints(MinEDt,distance,EminDT(index,:))
    
    drawnow
end
