clear
clc
%% MinCDE3D_height_wavelength
%   Generates 2 plots:
%   1. Max height v. time
%   2. Average wavelength v. time

%% Import
load('EminDT_redr.mat');
load('minDt_redr.mat');
load('time_redr.mat');
load('distance_redr.mat');

%% Process data
maxheightE = zeros(length(time_redr),1);
maxheightD = zeros(length(time_redr),1);

avglengthE = zeros(length(time_redr),1);
avglengthD = zeros(length(time_redr),1);

for i = 1:length(time_redr)
   pksE = findpeaks(EminDT_redr(i,:));
   maxheightE(i,1) = mean(pksE);
   localminE = findpeaks(-EminDT_redr(i,:));
   if length(pksE) == 1 %&& length(localminE) == 0
       avglengthE(i,1) = 0;
    elseif length(pksE) == 1
        [lowpksE lowpksEloc] = findpeaks(EminDT_redr(i,:));
        avgE = 0;
        for j = 1:length(lowpksEloc)-1
            avgE = avgE + (distance_redr(lowpksEloc(1,j+1),1)-distance_redr(lowpksEloc(1,j)));
        end
        avglengthE(i,1) = avgE/(length(lowpksEloc)-1);
    else
        [lowpksE lowpksEloc] = findpeaks(EminDT_redr(i,:));
        avgE = 0;
        for j = 1:length(lowpksEloc)-1
            avgE = avgE + (distance_redr(lowpksEloc(1,j+1),1)-distance_redr(lowpksEloc(1,j),1));
        end
        avglengthE(i,1) = avgE/(length(lowpksEloc)-1);
    end
    
    pksD = findpeaks(minDt_redr(i,:),'MinPeakDistance',4);
    maxheightD(i,1) = mean(pksD);
    if length(pksD) == 1
        avglengthD(i,1) = 0;
    else
        [lowpksD lowpksDloc] = findpeaks(minDt_redr(i,:),'MinPeakDistance',4);
        avgD = 0;
        for j = 1:length(lowpksDloc)-1
            avgD = avgD + (distance_redr(lowpksDloc(1,j+1),1)-distance_redr(lowpksDloc(1,j),1));
        end
        avglengthD(i,1) = avgD/(length(lowpksDloc)-1);
    end
end

%% Max height v. time plots
figure(3)
subplot(2,2,1)
plot(time_redr,maxheightE,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Height (N_{EminDt_{redr}}(t))','fontsize',14);
title('Average Maximum Height of EminDt for MinCDE 3D_{redr}','fontsize',16);

subplot(2,2,2)
plot(time_redr,maxheightD,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Height (N_{minDT_{redr}}(t))','fontsize',14);
title('Average Maximum Height of minDT for MinCDE 3D_{redr}','fontsize',16);

%% Average wavelength v. time
subplot(2,2,3)
plot(time_redr,avglengthE,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Average Wavelength (um))','fontsize',14);
title('Average Wavelength of EminDt for MinCDE 3D_{redr}','fontsize',16);

subplot(2,2,4)
plot(time_redr,avglengthD,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Average Wavelength (um)','fontsize',14);
title('Average Wavelength of minDT for MinCDE 3D_{redr}','fontsize',16);