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
ETotal = zeros(1,length(time));

for i = 2:length(time)
    [NumPks,PksInd] = findPeaks(EminDt(i,:),distance,Ethresh);
    ENumPks(i) = NumPks;
    EPksInd(i,1:length(PksInd)) = PksInd;
    EPksD(i, 1:length(PksInd)) = distance(PksInd);
    EPksDiff(i,1:length(PksInd)) = EPksD(i,1:length(PksInd))-[0 EPksD(i,1:length(PksInd)-1)];
    EPksHeight(i,1:length(PksInd)) = EminDt(i,PksInd);
    ESymVal(i) = symVal(EminDt(i,:),distance);
    ETotal(i) = rsum(EminDt(i,:)',distance');
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
DTotal = zeros(1,length(time));

for i = 2:length(time)
    [NumPks,PksInd] = findPeaks(minDt(i,:),distance,Dthresh);
    DNumPks(i) = NumPks;
    DPksInd(i,1:length(PksInd)) = PksInd;
    DPksD(i, 1:length(PksInd)) = distance(PksInd);
    DPksDiff(i,1:length(PksInd)) = DPksD(i,1:length(PksInd))-[0 DPksD(i,1:length(PksInd)-1)];
    DPksHeight(i,1:length(PksInd)) = minDt(i,PksInd);
    DSymVal(i) = symVal(minDt(i,:),distance);
    DTotal(i) = rsum(minDt(i,:)',distance');
end

DAvgHeight = sum(DPksHeight,2)'./DNumPks;
DAvgWL = sum(DPksDiff(:,2:floor(length(time)/4)),2)'./(DNumPks-1);
%% Period Calculations
species = {'EminDt','minDt'};
prop = {'PkHeight', 'Wavelength', 'SymVal', 'Total'};
% Period Averge over 1200 s
disp('Period Average Over 1200 s');
[epkm, epks] = findPeriod(EAvgHeight,time);
[ewlm, ewls] = findPeriod(EAvgWL, time);
[esymm, esyms] = findPeriod(ESymVal, time);
[etotm, etots] = findPeriod(ETotal, time);
[dpkm, dpks] = findPeriod(DAvgHeight,time);
[dwlm, dwls] = findPeriod(DAvgWL, time);
[dsymm, dsyms] = findPeriod(DSymVal, time);
[dtotm, dtots] = findPeriod(DTotal, time);
% Period Mean Calculations
disp('-----------------------------------------------------------');
disp('--Period Mean Over 1200s--');
pkheight = [epkm; dpkm];
wl = [ewlm; dwlm];
sym = [esymm; dsymm];
tot = [etotm; dtotm];
table(pkheight, wl, sym, tot, 'RowNames',species, 'VariableNames',prop)
% Period Std Calculations
disp('--Period Std Over 1200s--');
pkheights = [epks; dpks];
wls = [ewls; dwls];
syms = [esyms; dsyms];
tots = [etots; dtots];
table(pkheights, wls, syms, tots, 'RowNames',species, 'VariableNames',prop)

% Steady State Period Calculations
disp('Stead State Period Average Beyond 700 s');
[epkm, epks] = findPeriod(EAvgHeight,time,700);
[ewlm, ewls] = findPeriod(EAvgWL, time,700);
[esymm, esyms] = findPeriod(ESymVal, time,700);
[etotm, etots] = findPeriod(ETotal, time,700);
[dpkm, dpks] = findPeriod(DAvgHeight,time,700);
[dwlm, dwls] = findPeriod(DAvgWL, time,700);
[dsymm, dsyms] = findPeriod(DSymVal, time,700);
[dtotm, dtots] = findPeriod(DTotal, time,700);
% Period Mean Calculations
disp('-----------------------------------------------------------');
disp('--Period Mean Steady State--');
pkheight = [epkm; dpkm];
wl = [ewlm; dwlm];
sym = [esymm; dsymm];
tot = [etotm; dtotm];
table(pkheight, wl, sym, tot, 'RowNames',species, 'VariableNames',prop)
% Period Std Calculations
disp('--Period Std Steady State--');
pkheights = [epks; dpks];
wls = [ewls; dwls];
syms = [esyms; dsyms];
tots = [etots; dtots];
table(pkheights, wls, syms, tots, 'RowNames',species, 'VariableNames',prop)

%% Plot EminDt and minDt
%           Avg Peak Height Avg Wavelength Symmetry Value   Total Conc
%   EminDt          1             2              3              4
%   minDt           5             6              7              8

figure(1)
% EminDt - Average Peak Height over time
subplot(2,4,1)
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
xlabel('time(s)','fontsize',12);
ylabel('N_{EminDt}(t)','fontsize',12);
axis([0 1200 0 100]);
lgnd = legend(legendentries);
lgnd.FontSize = 8;
lgnd.Location = 'northwest'; 
lgnd.Orientation = 'vertical';

% EminDt - Average Wavelength Over Time
subplot(2,4,2)
EAvgWL(isnan(EAvgWL))=0;
plot(time, EAvgWL,'-','LineWidth',1);
title('Average Wavelength');
xlabel('time(s)','fontsize',12);
ylabel('distance (um)','fontsize',12);
axis([0 1200 0 2]);

% EminDt - Symmetry Indicator Over Time
subplot(2,4,3)
plot(time,ESymVal,'-','LineWidth',1);
title('Symmetry');
xlabel('time(s)','fontsize',12);
ylabel('Symmetry Number','fontsize',12);
axis([0 1200 -0.5 0.5]);

% EminDt - Chance in Total Number over Time
subplot(2,4,4)
plot(time,ETotal,'-','LineWidth',1);
title('Total Number of EMinDt');
xlabel('time(s)','fontsize',12);
ylabel('N_{EminDt}(t)','fontsize',12);
axis([0 1200 0 75])

% minDt - Average Peak Height
subplot(2,4,5)
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
xlabel('time(s)','fontsize',12);
ylabel('N_{minDt}(t)','fontsize',12);
axis([0 1200 0 100]);
lgnd = legend(legendentries);
lgnd.FontSize = 8;
lgnd.Location = 'northwest'; 
lgnd.Orientation = 'vertical';

% minDt - Average Wavelength
subplot(2,4,6)
DAvgWL(isnan(DAvgWL))=0;
plot(time, DAvgWL,'-','LineWidth',1);
title('Average Wavelength');
xlabel('time(s)','fontsize',12);
ylabel('distance (um)','fontsize',12);
axis([0 1200 0 2]);

% minDt - Symmetry Indicator
subplot(2,4,7)
plot(time,DSymVal,'-','LineWidth',1);
title('Symmetry');
xlabel('time(s)','fontsize',12);
ylabel('Symmetry Value');
axis([0 1200 -2 2]);

% minDt - Chance in Total Concentration over Time
subplot(2,4,8)
plot(time,DTotal,'-','LineWidth',1);
title('Total Number of MinDt');
xlabel('time(s)','fontsize',12);
ylabel('N_{minDt}(t)','fontsize',12);
axis([0 1200 0 100]);

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
%% findPeriod Function
%   Returns the expected value and standard deviation for the period 
%   Input:
%       data    must be a row vector
%       time    must be a column vector
%       startt  start time (double)
%       endt   end time (double)
function [permean, perstd] = findPeriod(data,time,startt,endt)
    switch nargin
        case 2
            startt = time(1);
            endt = time(length(time));
        case 3
            endt = time(length(time));
    end
    
    startind = find(time>= startt,1,'first');
    endind = find(time<= endt,1,'last');   
    
    [~, PksInd] = findPeaks(data(startind:endind),time(startind:endind));
    period = time(PksInd + startind - 1);
    perioddiff = period - [0; period(1:length(period)-1)];
    
    permean = mean(perioddiff(2:length(perioddiff)));
    perstd = std(perioddiff(2:length(perioddiff))); 
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
    symval = symval/length(data(center-range+1:center+range));
end
%% Reimann Sum Function
% Assumes data and time are column vectors
function [reisum] = rsum(data,time,startt,endt)
    % Check number of inputs.
    % Fill in unset optional values.
    switch nargin
        case 2
            startt = time(1);
            endt = time(length(time));
        case 3
            endt = time(length(time));
    end
    
    % Find time window; assign to vector called dist
    dist = time(find(time<=startt,1,'last'):find(time<=endt,1,'last'));
    
    % Calculate the difference in time and the total time covered
    diff = dist - [0; dist(1:size(dist,1)-1)];
    total = dist(size(dist,1)) - dist(1);
    
    % Preallocate Reimann Sum Component Vector
    reimann = zeros(size(dist,1),1);
    % Initialize last value
    reimann(length(reimann)) = data(length(reimann))*1/2*diff(length(reimann));
    
    for i = 1:length(reimann)-1
        reimann(i) = (1/2*diff(i)+1/2*diff(i+1))*data(i);
    end
    
    reisum = sum(reimann)/total;
end