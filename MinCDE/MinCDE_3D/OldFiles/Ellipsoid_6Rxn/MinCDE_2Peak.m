%% MinCDE_2Peak.m
%  Generates a plot of the average maximum height and wavelength for EminDt
%  and minDt
%% Initialize
clear
clc
%% Import
load('EminDt.mat');
load('minDt.mat');
load('time.mat');
load('distance.mat');

%% EminDt
% Define Peaks
peak1 = zeros(1,601);
peak1d = zeros(1,601);
peak2 = zeros(1,601);
peak2d = zeros(1,601);

for ctime = 2:length(time)
    data1 = EminDt(ctime,1:83);
    data2 = EminDt(ctime,84:167);
    maxval1 = max(data1);
    maxval2 = max(data2);
    index1 = find(data1==maxval1);
    index2 = find(data2==maxval2)+83;
    % Identifies whether there exist 1 or 2 peaks
    %   - abs(maxval1-maxval2) ensures there there exists a second peak in
    %   the second half of the cell of similar height (< 1 molecule apart)
    %   - index2-index1 ensures that the peaks are found at least 1 cell
    %   apart. 
    if(abs(maxval1-maxval2)<0.5 && (index2-index1) > 1)
        peak1(ctime) = maxval1;
        peak2(ctime) = maxval2;
        peak1d(ctime) = distance(index1);
        peak2d(ctime) = distance(index2);
    else
        peak1(ctime) = maxval1;
        peak1d(ctime) = distance(max(index1));
    end
end
% Assign peak names
Epk1 = peak1;
Epk2 = peak2;
Epk2(peak2 == 0) = NaN;
Epk1d = peak1d;
Epk2d = peak2d;

% Find average wavelength
Eavgwavel = zeros(1,length(time));
doublei = find(peak2d);
Eavgwavel(doublei) = peak2d(doublei)-peak1d(doublei);

% Find average peak height
Eavgpkh = peak1 + peak2;
Eavgpkh(doublei) = Eavgpkh(doublei)/2;

%% minDt
% Define Peaks
peak1 = zeros(1,601);
peak1d = zeros(1,601);
peak2 = zeros(1,601);
peak2d = zeros(1,601);

for ctime = 2:length(time)
    data1 = minDt(ctime,1:83);
    data2 = minDt(ctime,84:167);
    maxval1 = max(data1);
    maxval2 = max(data2);
    index1 = find(data1==maxval1);
    index2 = find(data2==maxval2)+83;
    % Identifies whether there exist 1 or 2 peaks
    %   - abs(maxval1-maxval2) ensures there there exists a second peak in
    %   the second half of the cell of similar height (< 1 molecule apart)
    %   - index2-index1 ensures that the peaks are found at least 1 cell
    %   apart. 
    if(abs(maxval1-maxval2) < 1 && index2-index1 > 1)
        peak1(ctime) = maxval1;
        peak2(ctime) = maxval2;
        peak1d(ctime) = distance(index1);
        peak2d(ctime) = distance(index2);
    else
        peak1(ctime) = maxval1;
        peak1d(ctime) = distance(max(index1));
    end
end
% Assign peak names
Dtpk1 = peak1;
Dtpk2 = peak2;
Dtpk2(peak2 == 0) = NaN;
Dtpk1d = peak1d;
Dtpk2d = peak2d;

% Find average wavelength
Dtavgwavel = zeros(1,length(time));
doublei = find(peak2d);
Dtavgwavel(doublei) = peak2d(doublei)-peak1d(doublei);

% Find average peak height
Dtavgpkh = peak1 + peak2;
Dtavgpkh(doublei) = Dtavgpkh(doublei)/2;

%% Figure 1 - Presence of 1 or 2 Peaks Over time
figure(1)
subplot(2,1,1)
hold on
plot(time,Epk1,'LineWidth',3);
plot(time,Epk2,'-','LineWidth',2);
xlabel('time(s)','fontsize',14);
ylabel('Peak Height (N_{EminDt}(t))','fontsize',14);
title('Peak Heights of EminDt','fontsize',16);

lgnd = legend('1 Peak','2 Peaks');
lgnd.FontSize = 11.5;
lgnd.Location = 'southeast'; 
lgnd.Orientation = 'vertical';

subplot(2,1,2)
hold on
plot(time,Dtpk1,'LineWidth',3);
plot(time,Dtpk2,'-','LineWidth',2);
xlabel('time(s)','fontsize',14);
ylabel('Peak Height (N_{minDt}(t))','fontsize',14);
title('Peak Heights of minDt','fontsize',16);

%% Figure 2 - Max height v. time plots
figure(2)
subplot(2,2,1)
plot(time,Eavgpkh,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Height (N_{EminDt}(t))','fontsize',14);
title('Average Maximum Height of EminDt for MinCDE','fontsize',16);

subplot(2,2,2)
plot(time,Dtavgpkh,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Height (N_{minDt}(t))','fontsize',14);
title('Average Maximum Height of minDT for MinCDE','fontsize',16);

%% Average wavelength v. time
subplot(2,2,3)
plot(time,Eavgwavel,'LineWidth',3);
xlabel('time (s)', 'fontsize',14);
ylabel('Average Wavelength (um))','fontsize',14);
title('Average Wavelength of EminDt for MinCDE 3D','fontsize',16);

subplot(2,2,4)
plot(time,Dtavgwavel,'LineWidth',3);
xlabel('time (s)','fontsize',14);
ylabel('Average Wavelength (um)','fontsize',14);
title('Average Wavelength of minDT for MinCDE 3D','fontsize',16);

%% Analysis of Temporal Changes in Average Height and Average Spatial Wavelength
% EminDt - Avg Wavelength
% Find index where zeros change to a positive value.
zerosi = find(Eavgwavel == 0);
patternf = zeros(0,length(time));
patternf(zerosi) = 1;
for r = 1:length(patternf)-1
    patternf(r) = patternf(r+1)-patternf(r);
end
zerosi=find(patternf~=0);
zerosi=zerosi(:,2:2:end);
zerosi = [1 zerosi];

% Find time at which the average height reaches a peak
%   Found by breaking down the Eavgwavel curve into segments separated by
%   periods of zeros
wl = zeros(1,length(zerosi)-1);
for ind = 1:length(zerosi)-1
    leftind = zerosi(ind);
    rightind = zerosi(ind+1);
    [maxv,index] = max(Eavgwavel(leftind:rightind));
    wl(ind)=time(index+leftind-1);
end

% Find average wavelength
avgwl = 0;
for a = 1:length(wl)-1
    avgwl = avgwl+wl(a+1)-wl(a);
end
Eavgwl = avgwl/(length(wl)-1)


% minDt - Avg Wavelength
% Find index where zeros change to a positive value.
zerosi = find(Dtavgwavel == 0);
patternf = zeros(0,length(time));
patternf(zerosi) = 1;
for r = 1:length(patternf)-1
    patternf(r) = patternf(r+1)-patternf(r);
end
zerosi=find(patternf~=0);
zerosi=zerosi(:,2:2:end);
zerosi = [1 zerosi];

% Find time at which the average height reaches a peak
%   Found by breaking down the Eavgwavel curve into segments separated by
%   periods of zeros
wl = zeros(1,length(zerosi)-1);
for ind = 1:length(zerosi)-1
    leftind = zerosi(ind);
    rightind = zerosi(ind+1);
    [maxv,index] = max(Dtavgwavel(leftind:rightind));
    wl(ind)=time(index+leftind-1);
end

% Find average wavelength
avgwl = 0;
for a = 1:length(wl)-1
    avgwl = avgwl+wl(a+1)-wl(a);
end
Dtavgwl = avgwl/(length(wl)-1)

