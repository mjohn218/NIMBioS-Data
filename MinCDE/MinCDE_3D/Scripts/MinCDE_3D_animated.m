%% MinCDE2D Animated Line

%% Import
load('minDt.mat');
load('EminDt.mat');
load('time.mat');
load('distance.mat');

%% Animation
figure(4)
clf
MinDt = animatedline('Color',[0 .75 1],'LineWidth',3);
MinEDt = animatedline('Color',[0 0 1],'LineWidth',3);

axis([0 6 0 200])
xlabel('Distance along long axis (um)','FontSize',14);
ylabel('N(x)','FontSize',14);
title('MinCDE (Ellipsoid, Original Parameters)','FontSize',16);
l=legend('MinDt','EminDT');
l.FontSize = 14;
for k = 1:length(time)*6
    clearpoints(MinDt)
    clearpoints(MinEDt)
    
    index = uint8(ceil(k/6));
    addpoints(MinDt,distance,minDt(index,:))
    addpoints(MinEDt,distance,EminDt(index,:))
    
    drawnow
end
