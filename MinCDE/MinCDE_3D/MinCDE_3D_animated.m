%% MinCDE3D Animated Line

%% Import
load('minDt.mat');
load('time.mat');
load('distance.mat');

%% Animation
figure(4)
clf
MinCDE3D = animatedline;
axis([0 6 0 100])
xlabel('Distance along long axis (um)','FontSize',14);
ylabel('N_{MinDt}(x)','FontSize',14);
title('Changes in MinDt over Time in MinCDE3D','FontSize',16);
MinCDE3D.LineWidth = 3;
MinCDE3D.Color = [0 .75 1];
for k = 1:length(time)*6
    clearpoints(MinCDE3D)
    index = uint8(ceil(k/6));
    addpoints(MinCDE3D,distance,minDt(index,:))
    drawnow
end