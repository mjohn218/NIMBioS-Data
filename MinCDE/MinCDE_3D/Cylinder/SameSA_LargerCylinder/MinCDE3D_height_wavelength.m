clear
clc
%% MinCDE3D_height_wavelength
%   Generates 2 plots:
%   1. Max height v. time
%   2. Average wavelength v. time

%% Import
load('EminDt_enlarged.mat');
load('minDt_enlarged.mat');
load('time_enlarged.mat');
load('distance_enlarged.mat');

%% Process data
maxheightE = zeros(length(time_enlarged),1);
maxheightD = zeros(length(time_enlarged),1);

avglengthE = zeros(length(time_enlarged),1);
avglengthD = zeros(length(time_enlarged),1);

for i = 1:length(time_enlarged)
   pksE = findpeaks(EminDT_enlarged(i,:));
   maxheightE(i,1) = mean(pksE);
   localminE = findpeaks(-EminDT_enlarged(i,:));
   if length(pksE) == 1 %&& length(localminE) == 0
       avglengthE(i,1) = 0;
    elseif length(pksE) == 1
        [lowpksE lowpksEloc] = findpeaks(EminDT_enlarged(i,:));
        avgE = 0;
        for j = 1:length(lowpksEloc)-1
            avgE = avgE + (distance_enlarged(lowpksEloc(1,j+1),1)-distance_enlarged(lowpksEloc(1,j)));
        end
        avglengthE(i,1) = avgE/(length(lowpksEloc)-1);
    else
        [lowpksE lowpksEloc] = findpeaks(EminDT_enlarged(i,:));
        avgE = 0;
        for j = 1:length(lowpksEloc)-1
            avgE = avgE + (distance_enlarged(lowpksEloc(1,j+1),1)-distance_enlarged(lowpksEloc(1,j),1));
        end
        avglengthE(i,1) = avgE/(length(lowpksEloc)-1);
    end
    
    pksD = findpeaks(minDt_enlarged(i,:),'MinPeakDistance',4);
    maxheightD(i,1) = mean(pksD);
    if length(pksD) == 1
        avglengthD(i,1) = 0;
    else
        [lowpksD lowpksDloc] = findpeaks(minDt_enlarged(i,:),'MinPeakDistance',4);
        avgD = 0;
        for j = 1:length(lowpksDloc)-1
            avgD = avgD + (distance_enlarged(lowpksDloc(1,j+1),1)-distance_enlarged(lowpksDloc(1,j),1));
        end
        avglengthD(i,1) = avgD/(length(lowpksDloc)-1);
    end
end

%% Max height v. time plots
figure(3)
subplot(2,2,1)
plot(time_enlarged,maxheightE,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Height (N_{EminDt}(t))','fontsize',14);
title('Average Maximum Height of EminDt for MinCDE 3D - Enlarged','fontsize',16);

subplot(2,2,2)
plot(time_enlarged,maxheightD,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Height (N_{minDT}(t))','fontsize',14);
title('Average Maximum Height of minDT for MinCDE 3D - Enlarged','fontsize',16);

%% Average wavelength v. time
subplot(2,2,3)
plot(time_enlarged,avglengthE,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Average Wavelength (um))','fontsize',14);
title('Average Wavelength of EminDt for MinCDE 3D - Enlarged','fontsize',16);

subplot(2,2,4)
plot(time_enlarged,avglengthD,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Average Wavelength (um)','fontsize',14);
title('Average Wavelength of minDT for MinCDE 3D-Enlarged','fontsize',16);