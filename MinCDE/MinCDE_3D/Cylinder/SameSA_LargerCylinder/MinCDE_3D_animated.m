%% MinCDE2D Animated Line

%% Import
load('minDt_enlarged.mat');
load('EminDT_enlarged.mat');
load('time_enlarged.mat');
load('distance_enlarged.mat');

%% Animation
figure(4)
clf
MinDt = animatedline('Color',[0 .75 1],'LineWidth',3);
MinEDt = animatedline('Color',[0 0 1],'LineWidth',3);

axis([0 2.6 0 200])
xlabel('Distance along long axis (um)','FontSize',14);
ylabel('N(x)','FontSize',14);
title('Changes in Number of Particles over Time in MinCDE3D','FontSize',16);
l=legend('MinDt','EminDT');
l.FontSize = 14;
for k = 1:length(time_enlarged)*8
    clearpoints(MinDt)
    clearpoints(MinEDt)
    
    index = uint8(ceil(k/8));
    addpoints(MinDt,distance_enlarged,minDt_enlarged(index,:))
    addpoints(MinEDt,distance_enlarged,EminDT_enlarged(index,:))
    
    drawnow
end
