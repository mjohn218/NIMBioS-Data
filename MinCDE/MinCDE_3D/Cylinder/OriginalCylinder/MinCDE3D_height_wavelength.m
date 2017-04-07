clear
clc
%% MinCDE3D_height_wavelength
%   Generates 2 plots:
%   1. Max height v. time
%   2. Average wavelength v. time

%% Import
load('EminDt.mat');
load('minDt.mat');
load('time.mat');
load('distance.mat');

%% Process data
maxheightE = zeros(length(time),1);
maxheightD = zeros(length(time),1);

avglengthE = zeros(length(time),1);
avglengthD = zeros(length(time),1);

for i = 1:length(time)
   pksE = findpeaks(EminDt(i,:));
   maxheightE(i,1) = mean(pksE);
   localminE = findpeaks(-EminDt(i,:));
   if length(pksE) == 1 %&& length(localminE) == 0
       avglengthE(i,1) = 0;
    elseif length(pksE) == 1
        [lowpksE lowpksEloc] = findpeaks(EminDt(i,:));
        avgE = 0;
        for j = 1:length(lowpksEloc)-1
            avgE = avgE + (distance(lowpksEloc(1,j+1),1)-distance(lowpksEloc(1,j)));
        end
        avglengthE(i,1) = avgE/(length(lowpksEloc)-1);
    else
        [lowpksE lowpksEloc] = findpeaks(EminDt(i,:));
        avgE = 0;
        for j = 1:length(lowpksEloc)-1
            avgE = avgE + (distance(lowpksEloc(1,j+1),1)-distance(lowpksEloc(1,j),1));
        end
        avglengthE(i,1) = avgE/(length(lowpksEloc)-1);
    end
    
    pksD = findpeaks(minDt(i,:),'MinPeakDistance',1);
    maxheightD(i,1) = mean(pksD);
    if length(pksD) == 1
        avglengthD(i,1) = 0;
    else
        [lowpksD lowpksDloc] = findpeaks(minDt(i,:),'MinPeakDistance',1);
        avgD = 0;
        for j = 1:length(lowpksDloc)-1
            avgD = avgD + (distance(lowpksDloc(1,j+1),1)-distance(lowpksDloc(1,j),1));
        end
        avglengthD(i,1) = avgD/(length(lowpksDloc)-1);
    end
end

maxheightD(isnan(maxheightD)) = 0;
maxheightE(isnan(maxheightE)) = 0;
avglengthD(isnan(avglengthD)) = 0;
avglengthE(isnan(avglengthE)) = 0;

%% Max height v. time plots
figure(3)
subplot(2,2,1)
plot(time,maxheightE,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Height (N_{EminDt}(t))','fontsize',14);
title('Average Maximum Height of EminDt for MinCDE 3D','fontsize',16);

subplot(2,2,2)
plot(time,maxheightD,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Height (N_{minDT}(t))','fontsize',14);
title('Average Maximum Height of minDT for MinCDE 3D','fontsize',16);

%% Average wavelength v. time
subplot(2,2,3)
plot(time,avglengthE,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Average Wavelength (um))','fontsize',14);
title('Average Wavelength of EminDt for MinCDE 3D','fontsize',16);

subplot(2,2,4)
plot(time,avglengthD,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Average Wavelength (um)','fontsize',14);
title('Average Wavelength of minDT for MinCDE 3D','fontsize',16);