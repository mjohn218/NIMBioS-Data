clear
clc
%% MinCDE3D_height_wavelength
%   Generates 2 plots:
%   1. Max height v. time
%   2. Average wavelength v. time

%% Import
load('EminDt_R26.mat');
load('minDt_R26.mat');
load('time_R26.mat');
load('distance_R26.mat');

%% Process data
maxheightE = zeros(length(time_R26),1);
maxheightD = zeros(length(time_R26),1);

avglengthE = zeros(length(time_R26),1);
avglengthD = zeros(length(time_R26),1);

for i = 1:length(time_R26)
   pksE = findpeaks(EminDt_R26(i,:),'MinPeakDistance',4);
   maxheightE(i,1) = mean(pksE);
   localminE = findpeaks(-EminDt_R26(i,:),'MinPeakDistance',4);
   if length(pksE) == 1 %&& length(localminE) == 0
       avglengthE(i,1) = 0;
    elseif length(pksE) == 1
        [lowpksE lowpksEloc] = findpeaks(EminDt_R26(i,:),'MinPeakDistance',4);
        avgE = 0;
        for j = 1:length(lowpksEloc)-1
            avgE = avgE + (distance_R26(lowpksEloc(1,j+1),1)-distance_R26(lowpksEloc(1,j)));
        end
        avglengthE(i,1) = avgE/(length(lowpksEloc)-1);
    else
        [lowpksE lowpksEloc] = findpeaks(EminDt_R26(i,:));
        avgE = 0;
        for j = 1:length(lowpksEloc)-1
            avgE = avgE + (distance_R26(lowpksEloc(1,j+1),1)-distance_R26(lowpksEloc(1,j),1));
        end
        avglengthE(i,1) = avgE/(length(lowpksEloc)-1);
    end
    
    pksD = findpeaks(minDt_R26(i,:),'MinPeakDistance',4);
    maxheightD(i,1) = mean(pksD);
    if length(pksD) == 1
        avglengthD(i,1) = 0;
    else
        [lowpksD lowpksDloc] = findpeaks(minDt_R26(i,:),'MinPeakDistance',4);
        avgD = 0;
        for j = 1:length(lowpksDloc)-1
            avgD = avgD + (distance_R26(lowpksDloc(1,j+1),1)-distance_R26(lowpksDloc(1,j),1));
        end
        avglengthD(i,1) = avgD/(length(lowpksDloc)-1);
    end
end

%% Max height v. time plots
figure(3)
subplot(2,2,1)
plot(time_R26,maxheightE,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Height (N_{EminDt_{R26}}(t))','fontsize',14);
title('Average Maximum Height of EminDt for MinCDE 3D_{R26}','fontsize',16);

subplot(2,2,2)
plot(time_R26,maxheightD,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Height (N_{minDT_{R26}}(t))','fontsize',14);
title('Average Maximum Height of minDT for MinCDE 3D_{R26}','fontsize',16);

%% Average wavelength v. time
subplot(2,2,3)
plot(time_R26,avglengthE,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Average Wavelength (um))','fontsize',14);
title('Average Wavelength of EminDt for MinCDE 3D_{R26}','fontsize',16);

subplot(2,2,4)
plot(time_R26,avglengthD,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Average Wavelength (um)','fontsize',14);
title('Average Wavelength of minDT for MinCDE 3D_{R26}','fontsize',16);