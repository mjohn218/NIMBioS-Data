%Read in a concatenated file where it reports A then R in each column.
%So min size is 2 columns for AR. 
function[sepsA, sepsR, lags, locA, locR]=calc_peak_seps(time, AR)

%This is set up to ignore peaks that are too narrow (MinPeakWidth)
%It is also set up based on what we know: the peaks are ~20s apart. That
%way it will ignore any small peaks that occur.
%Or can use MinPeakHeight

[len, N]=size(AR);
Ntraj=N/2;
sepsA=[];
sepsR=[];
lags=[];
for i=1:1:Ntraj
    
    ind=(i-1)*2+1
    [pksA,locA]=findpeaks(AR(:,ind),time, 'MinPeakWidth',2,'MinPeakHeight',500);
    [pksR,locR]=findpeaks(AR(:,ind+1),time, 'MinPeakWidth',2,'MinPeakHeight',500);
    size(pksA)
    size(pksR)
    %mA(i)=mean(diff(locA));
    %stdA(i)=std(diff(locA));
    sepsA=[sepsA;diff(locA)];
   % medianA(i)=median(diff(locA));
    sepsR=[sepsR;diff(locR)];
    if(length(locR)==length(locA))
        lags=[lags;(locR-locA)];
    else
        display('Num peaks mismatch')
        lags=[lags;(locR-locA(1:end-1))];
    end
end

        