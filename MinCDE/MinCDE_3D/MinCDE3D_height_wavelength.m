clear
clc
%% MinCDE3D_height_wavelength
%   Generates 2 plots:
%   1. Max height v. time
%   2. Average wavelength v. time

%% Import
%load('EminDT.mat');
load('minDt.mat');
load('time.mat');
load('distance.mat');

%% Process data
%maxheightE = zeros(length(time),1);
maxheightD = zeros(length(time),1);

%avglengthE = zeros(length(time),1);
avglengthD = zeros(length(time),1);

for i = 1:length(time)
%    pksE = findpeaks(EminDT(i,:));
%    maxheightE(i,1) = mean(pksE);
%    localminE = findpeaks(-EminDT(i,:));
%    if length(pksE) == 1 %&& length(localminE) == 0
%        avglengthE(i,1) = 0;
%     elseif length(pksE) == 1
%         [lowpksE lowpksEloc] = findpeaks(EminDT(i,:));
%         avgE = 0;
%         for j = 1:length(lowpksEloc)-1
%             avgE = avgE + (distance(lowpksEloc(1,j+1),1)-distance(lowpksEloc(1,j)));
%         end
%         avglengthE(i,1) = avgE/(length(lowpksEloc)-1);
%     else
%         [lowpksE lowpksEloc] = findpeaks(EminDT(i,:));
%         avgE = 0;
%         for j = 1:length(lowpksEloc)-1
%             avgE = avgE + (distance(lowpksEloc(1,j+1),1)-distance(lowpksEloc(1,j),1));
%         end
%         avglengthE(i,1) = avgE/(length(lowpksEloc)-1);
%     end
    
    pksD = findpeaks(minDt(i,:));
    maxheightD(i,1) = mean(pksD);
    if length(pksD) == 1
        avglengthD(i,1) = 0;
    else
        [lowpksD lowpksDloc] = findpeaks(minDt(i,:));
        avgD = 0;
        for j = 1:length(lowpksDloc)-1
            avgD = avgD + (distance(lowpksDloc(1,j+1),1)-distance(lowpksDloc(1,j),1));
        end
        avglengthD(i,1) = avgD/(length(lowpksDloc)-1);
    end
end

%% Max height v. time plots
figure(3)
% subplot(2,2,1)
% plot(time,maxheightE,'LineWidth',3);
% xlabel('time (s)', 'fontsize',14);
% ylabel('Height (N_{EminDT}(t))','fontsize',14);
% title('Average Maximum Height of EminDT for MinCDE 3D','fontsize',16);

subplot(2,1,1)
plot(time,maxheightD,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Height (N_{minDT}(t))','fontsize',14);
title('Average Maximum Height of minDT for MinCDE 3D','fontsize',16);

%% Average wavelength v. time
% subplot(2,1,2)
% plot(time,avglengthE,'LineWidth',3);
% xlabel('time (s)', 'fontsize',14);
% ylabel('Average Wavelength (um))','fontsize',14);
% title('Average Wavelength of EminDT for MinCDE 3D','fontsize',16);

subplot(2,1,2)
plot(time,avglengthD,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Average Wavelength (um)','fontsize',14);
title('Average Wavelength of minDT for MinCDE 3D','fontsize',16);