%% MinCDE_Peaks

%% Initialize
clear
clc
close all

%% Import
load('distance.mat');
load('time.mat');
load('EminDt.mat');
load('minDt.mat');

%% EminDt
ENumPks = zeros(1,length(time));
EPksInd = zeros(length(time),floor(length(time)/4));
EPksD = zeros(length(time),floor(length(time)/4));
EPksDiff = zeros(length(time),floor(length(time)/4));
EPksHeight = zeros(length(time),floor(length(time)/4));
ESymVal = zeros(1,length(time));
Ethresh = max(max(EminDt))*2.5e-4;

for i = 2:length(time)
    [NumPks,PksInd] = findPeaks(EminDt(i,:),distance,Ethresh);
    ENumPks(i) = NumPks;
    EPksInd(i,1:length(PksInd)) = PksInd;
    EPksD(i, 1:length(PksInd)) = distance(PksInd);
    EPksDiff(i,1:length(PksInd)) = EPksD(i,1:length(PksInd))-[0 EPksD(i,1:length(PksInd)-1)];
    EPksHeight(i,1:length(PksInd)) = EminDt(i,PksInd);
    ESymVal(i) = symVal(EminDt(i,:),distance);
end

EAvgHeight = sum(EPksHeight,2)'./ENumPks;
EAvgWL = sum(EPksDiff(:,2:floor(length(time)/4)),2)'./(ENumPks-1);
%% minDt
DNumPks = zeros(1,length(time));
DPksInd = zeros(length(time),floor(length(time)/4));
DPksD = zeros(length(time),floor(length(time)/4));
DPksDiff = zeros(length(time),floor(length(time)/4));
DPksHeight = zeros(length(time),floor(length(time)/4));
DSymVal = zeros(1,length(time));
Dthresh = max(max(minDt))*2.5e-4;

for i = 2:length(time)
    [NumPks,PksInd] = findPeaks(minDt(i,:),distance,Dthresh);
    DNumPks(i) = NumPks;
    DPksInd(i,1:length(PksInd)) = PksInd;
    DPksD(i, 1:length(PksInd)) = distance(PksInd);
    DPksDiff(i,1:length(PksInd)) = DPksD(i,1:length(PksInd))-[0 DPksD(i,1:length(PksInd)-1)];
    DPksHeight(i,1:length(PksInd)) = minDt(i,PksInd);
    DSymVal(i) = symVal(minDt(i,:),distance);
end

DAvgHeight = sum(DPksHeight,2)'./DNumPks;
DAvgWL = sum(DPksDiff(:,2:floor(length(time)/4)),2)'./(DNumPks-1);

%% Plot EminDt and minDt
%           Avg Peak Height    Avg Wavelength      Symmetry Value
%   EminDt          1                 2                  3
%   minDt           4                 5                  6

figure(1)
% EminDt - Average Peak Height over time
subplot(2,3,1)
hold on
legendentries = cell(1,1);
% Plot the first curve
plot(time, EAvgHeight,'-','LineWidth',1);
legendentries{1,1} = '1 peak';
for p = 2:max(ENumPks)
    if(length(find(ENumPks==p))~=0)
        data = EAvgHeight;
        data(ENumPks ~= p) = NaN;
        plot(time,data,'-','LineWidth',1);
        legendentries = [legendentries strcat(num2str(p),' peaks')];
    end
end
title('Average Peak Height');
ylabel('N_{EminDt}(t)','fontsize',12);
axis([0 1200 0 1500]);
lgnd = legend(legendentries);
lgnd.FontSize = 8;
lgnd.Location = 'northwest'; 
lgnd.Orientation = 'vertical';
% Period calculation
[PNumPks, PPksInd] = findPeaks(EAvgHeight, time);
period = time(PPksInd);
perioddiff = period - [0; period(1:length(period)-1)];
period = sum(perioddiff)/length(perioddiff);
annotation('textbox',[.27 .825 .1 .1],'String',strcat('Period: ', num2str(period),' s'),'LineStyle','none');


% EminDt - Average Wavelength Over Time
subplot(2,3,2)
EAvgWL(isnan(EAvgWL))=0;
plot(time, EAvgWL,'-','LineWidth',1);
title('Average Wavelength');
xlabel('time (s)','fontsize',12);
ylabel('distance (um)','fontsize',12);
axis([0 1200 0 6]);
% Period calculation
[PNumPks, PPksInd] = findPeaks(EAvgWL, time);
period = time(PPksInd);
perioddiff = period - [0; period(1:length(period)-1)];
period = sum(perioddiff)/length(perioddiff);
annotation('textbox',[.55 .825 .1 .1],'String',strcat('Period: ', num2str(period),' s'),'LineStyle','none');


% EminDt - Symmetry Indicator Over Time
subplot(2,3,3)
plot(time,ESymVal,'-','LineWidth',1);
title('Symmetry');
ylabel('Symmetry Number','fontsize',12);
axis([0 1200 -3e4 3e4]);
% Period calculation
[PNumPks, PPksInd] = findPeaks(ESymVal, time);
period = time(PPksInd);
perioddiff = period - [0; period(1:length(period)-1)];
period = sum(perioddiff)/length(perioddiff);
annotation('textbox',[.83 .825 .1 .1],'String',strcat('Period: ', num2str(period),' s'),'LineStyle','none');



% minDt - Average Peak Height
subplot(2,3,4)
hold on
legendentries = cell(1,max(DNumPks));
% Plot the first curve
plot(time, DAvgHeight,'-','LineWidth',1);
legendentries{1,1} = '1 peak';
for p = 2:max(DNumPks)
    data = DAvgHeight;
    data(DNumPks ~= p) = NaN;
    plot(time,data,'-','LineWidth',1);
    legendentries{1,p} = strcat(num2str(p),' peaks');
end
title('Average Peak Height');
ylabel('N_{minDt}(t)','fontsize',12);
axis([0 1200 0 400]);
lgnd = legend(legendentries);
lgnd.FontSize = 8;
lgnd.Location = 'northwest'; 
lgnd.Orientation = 'vertical';
% Period calculation
[PNumPks, PPksInd] = findPeaks(DAvgHeight, time);
period = time(PPksInd);
perioddiff = period - [0; period(1:length(period)-1)];
period = sum(perioddiff)/length(perioddiff);
annotation('textbox',[.27 .35 .1 .1],'String',strcat('Period: ', num2str(period),' s'),'LineStyle','none');


% minDt - Average Wavelength
subplot(2,3,5)
DAvgWL(isnan(DAvgWL))=0;
plot(time, DAvgWL,'-','LineWidth',1);
title('Average Wavelength');
xlabel('time(s)','fontsize',12);
axis([0 1200 0 6]);
[PNumPks, PPksInd] = findPeaks(DAvgWL, time);
period = time(PPksInd);
perioddiff = period - [0; period(1:length(period)-1)];
period = sum(perioddiff)/length(perioddiff);
annotation('textbox',[.55 .35 .1 .1],'String',strcat('Period: ', num2str(period),' s'),'LineStyle','none');


% minDt - Symmetry Indicator
subplot(2,3,6)
plot(time,DSymVal,'-','LineWidth',1);
title('Symmetry');
ylabel('Symmetry Value');
axis([0 1200 -5e-3 5e-3]);
[PNumPks, PPksInd] = findPeaks(DSymVal, time);
period = time(PPksInd);
perioddiff = period - [0; period(1:length(period)-1)];
period = sum(perioddiff)/length(perioddiff);
annotation('textbox',[.83 .35 .1 .1],'String','Period: NaN','LineStyle','none');


%% FindPeaksFunction
%   Takes the concentration data and cell distance vector
%   Returns the number of peaks and the index of each peak.
function [NumPks,PksInd] = findPeaks(data, distance,threshold)
    % Check number of inputs.
    % Fill in unset optional values.
    switch nargin
        case 2
            threshold = max(data)*2.5e-4;
    end
    % Look at the difference between the two points
    diff = data - [0 data(1:length(data)-1)];
    
    % If the difference is significant (> threshold), assign +1 for positive 
    % differences and -1 for negative differences.
    % Otherwise, take on the previous trend value.
    if(diff(1)>0)
        diff(1)=1;
    else
        diff(1)=-1;
    end
    
    for i = 2:length(diff)
        if(diff(i) > threshold)
            diff(i)=1;
        elseif(diff(i)<-threshold)
            diff(i)=-1;
        else
            diff(i)=diff(i-1);
        end
    end
    
    % Identify points where the trend switches from +1 to -1.
    pks = zeros(1,length(distance));
    for i = 1:length(distance)-1
        if(diff(i) > 0 && diff(i+1)<0)
            pks(i)=1;
        end
    end

    % Return number of peaks and their indices.
    PksInd = find(pks);
    NumPks = length(PksInd);
end

%% Symmetry Measure
%   Returns a value indicating the symmetry of the cell at a given
%   timepoint. 
%   Calculated by taking the left values and subtracting the values on the
%   right. A negative value indicates that the cell is asymmetric with higher
%   values on the right side of the cell, and a positive value vice versa.
function symval = symVal(data,distance)
    maxd = max(distance);
    center = find(distance<= maxd/2,1,'last');
    range = min(center,length(distance)-center);
    symval = sum(data(center-range+1:center)-fliplr(data(center+1:center+range)));
end