clear
clc
close all
%% MinCDE2D_height_wavelength
%   Generates 2 plots:
%   1. Max height v. time
%   2. Average wavelength v. time

%% Import
load('EminDT.mat');
load('minDt.mat');
load('time.mat');
load('distance.mat');

%% Process data
maxheightE = zeros(length(time),1);
maxheightD = zeros(length(time),1);

avglengthE = zeros(length(time),1);
avglengthD = zeros(length(time),1);

for i = 1:length(time)

    wlt(i) = findpeaksfft(EminDT(i,:),distance);
       
end
figure
plot(time,wlt,'LineWidth',3);

% %% Max height v. time plots
% figure(3)
% subplot(2,2,1)
% plot(time,maxheightE,'LineWidth',3);
% xlabel('time (s)', 'fontsize',14);
% ylabel('Height (N_{EminDT}(t))','fontsize',14);
% title('Average Maximum Height of EminDT for MinCDE 2D','fontsize',16);
% 
% subplot(2,2,2)
% set(gca,'fontsize',18, 'fontweight','bold');
% plot(time,maxheightD,'LineWidth',3);
% xlabel('time (s)','fontsize',14);
% ylabel('Height (N_{minDT}(t))','fontsize',14);
% title('Average Maximum Height of minDT for MinCDE 2D','fontsize',16);
% 
% %% Average wavelength v. time
% subplot(2,2,3)
% plot(time,avglengthE,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Average Wavelength (um)','fontsize',14);
title('Average Wavelength of EminDT for MinCDE 2D','fontsize',16);
% 
% subplot(2,2,4)
% set(gca,'fontsize',18, 'fontweight','bold');
% plot(time,avglengthD,'LineWidth',3);
% xlabel('time (s)','fontsize',14);
% ylabel('Average Wavelength (um)','fontsize',14);
% title('Average Wavelength of minDT for MinCDE 2D','fontsize',16);