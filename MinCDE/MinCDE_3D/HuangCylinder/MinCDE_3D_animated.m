%% MinCDE2D Animated Line

%%


%% Import
load('minDt.mat');
load('EminDt.mat');
load('time.mat');
load('distance.mat');

%% Animation
figure(7)
clf
MinDt = animatedline('Color',[0 .75 1],'LineWidth',3);
MinEDt = animatedline('Color',[0 0 1],'LineWidth',3);

axis([0 4 0 500])
xlabel('Distance along long axis (um)','FontSize',14);
ylabel('N(x) (molecules/um^2)','FontSize',14);
title('MinCDE (Cylinder, long edge only)','FontSize',16);
l=legend('MinDt','EminDT');
l.FontSize = 14;
%Frame repeat here is 6
frame_rep=2;
v=VideoWriter('matlab_movie.avi');
open(v);
for k = 1:length(time)*frame_rep
    clearpoints(MinDt)
    clearpoints(MinEDt)
    
    index = uint8(ceil(k/frame_rep));
    addpoints(MinDt,distance,minDt(index,:))
    addpoints(MinEDt,distance,EminDt(index,:))
    
    drawnow
    M(k)=getframe(gcf);
    writeVideo(v, M(k));
end
close(v);