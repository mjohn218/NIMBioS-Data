%takes in a file with first col time (here in us! it is multiplied by 10^-6 to convert to s), next cols all copy
%numbers. Averages over all trajectories. 
%Important: For fitting, specify Volume and N0 inside of fitcrowd and
%fitcrowdlog.m!!!
%whichline=1, new plot. whichline>1, add a line to existing plot, color
%chosen according to value of whichline (<62).

%output: StoreValues has fits to A=A0exp(-B0t*kon), using the exp function
%and fitting the logy to a line (loses all y=0 data!!!).
%kon is initial rate, in units of nm^3/s.
%rawDataAvg is [time (s), mean(Copy)_over_alltraj, std(Copy)_over_all_traj


function[StoreValues, strInterval,rawDataAvg, beta, betalog]=analyze_crowding_commented(fprefix,fend, nsets, kon, figNum, whichline)

%Read in Input files
  %fname=sprintf('continuum_K%g_dt0.1_nBndPair.dat',Kd);
%fname=sprintf('membind_%g_set%d_nA.dat',Kd, setNum);
%box=300;k
nminLines=1E8;
nmaxLines=0;
%FIrst get the length of shortest file.
totcol=0;
for i=1:1:nsets
    fname=sprintf('%s_%s',fprefix,  fend);
    data=load(fname);
    [nlines,ncol]=size(data);
    totcol=totcol+ncol-1;
    if(nlines<nminLines)
        nminLines=nlines
    end
    if(nlines>nmaxLines)
        nmaxLines=nlines
    end
end
startc=1;
display('longest file')
nmaxLines
display('shortest file')
nminLines
time=data(:,1);
[a,b]=find(time<0);
lstart=1;
if(length(a)>0)
    display('skip negative times')
    
    lstart=a(end)+1;
    lstart
end
ltot=nminLines-lstart+1;
nBoundPair=zeros(ltot, totcol);
for i=1:1:nsets
    fname=sprintf('%s_set%g_%s',fprefix, i, fend);
    data=load(fname);
    [nlines,ncol]=size(data);
    endc=ncol-2+startc
    time=data(lstart:nminLines, 1);
   nBoundPair(:, startc:endc)=data(lstart:nminLines, 2:end);
    startc=endc+1
end



readCnotA=0;
loopCol=1;%use 2 if nloops is present!!

if(nBoundPair(1,2:end)==0)
    readCnotA=1;
    display('READING IN C (BOUND STATE) VALUES');
end
if(nBoundPair(:,3)==0)
    loopCol=2;
    display('WARNING: EVERY OTHER COLUMN IS SKIPPED, ASSUMES NLOOP COL PRESENT');
end
display('size')
size(nBoundPair)
%fname2=sprintf('memCltDrot_K%g_dt%g_meanComplexSize.dat',Kd, deltat);
%meanComplexSize=load(fname2);
%[r1, c1]=size(meanComplexSize);
%for i=2:1:c1
%    ind=find(isnan(nBoundPair(:,i)));
%    meanComplexSize(ind,i)=zeros(length(ind),1);
%    clear ind;
%end

%Mean over all 20 trajectories of meanComplexSize
%mtot=[meanComplexSize(:,1)*1E-6, mean(meanComplexSize(:,2:end)')', std(meanComplexSize(:,2:end)')'];
 %time=nBoundPair(:,1)*1E-6; %move time from us to s
 time=time*1E-6;
TimeLength=time(end);
tstart=TimeLength/5;
tfinish=tstart*2.5;
%Calculate m4eans from 2-5s and from 5s to the end.
%tmp=find(mtot(:,1)>tstart);
%%tmp=find(mtot(:,1)<tfinish);
%i2=tmp(end);
%ComplexSize2to5=mean(mtot(i1:i2,2));
%ComplexSize5toend=mean(mtot(i2:end,2));
%stdCSize2to5=mean(mtot(i1:i2,3));
%stdCSize5toend=mean(mtot(i2:end,3));

%Keq=kon/koff;%. 
%Calculate theoretical value of Ceq, assuming A+A<->C

%Use average ffrom the simulation itself
%at=sum(nBoundPair(:,col)')';
%ct=(at(1)-at)/2;

%For nbound pairs, it contains iters
[r, cols]=size(nBoundPair);

if(readCnotA==1)
    nBoundPair=[100-nBoundPair];
end

	ct=mean(nBoundPair(:,2:loopCol:end)')'; %average Bound legs over all trajectories
	display('final ct value');
	ct(end)
	complext=ct(1)-ct;
		stdct=std(nBoundPair(:,2:loopCol:end)')';
        rawDataAvg=[time, ct, stdct];

%Calculate Simulated value of Ceq, given A+A<->C but with multi-domain
%trimers (some legs cannot bind due to structural constraints). 



%Number of closed loops formed (hexagons). 
%nloopmean=mean(nBoundPair(:,3:2:end)')';%Average N loops
%nloopstd=std(nBoundPair(:,3:2:end)')';
if(loopCol==1)
    ntraj=(cols);
else
    ntraj=(cols-1)/2
end
%newt=linspace(time(1), time(end), 1E5);
%delt=newt(2)-newt(1);


%PLOT TIME DEPENDENCE
colormat=[0 1 1; 0 0 1; 0.494117647409439 0.184313729405403 0.556862771511078;...
    1 0 1; 0.466666668653488 0.674509823322296 0.18823529779911;...
    1 0 0; 1 0.7 0; 0.5 0.5 0.5; 0 0 0];

f1=figure(figNum)
%Reduce the number of data points plotted for large arrays.
nskip=5;
display('ntraj')
beta=zeros(ntraj,1);
betalog=zeros(ntraj,1);
    for i=1:ntraj
      
        [beta(i),R]=nlinfit(time,nBoundPair(:,i),@fitcrowd, kon);
     %   [betalog(i),R]=nlinfit(time,log(nBoundPair(:,i)),@fitcrowdlog, kon);
    end

%Allow for solid (ls=1) or dashed (ls=2) linestyles.
lineStyle{1}='-';
lineStyle{2}='--';
	     ls=1;
         colormat=colormap;
if(whichline==1)
    ax1=axes('Parent',f1);
 hold(ax1)


    semilogy(time(1:nskip:end), ct(1:nskip:end),lineStyle{ls},'LineWidth',5,'Color','r');
    % Create xlabel
    
xlabel('time (s)','FontWeight','bold');
    % Create ylabel
ylabel('A(t)','FontWeight','bold');
set(ax1, 'LineWidth',3,'FontSize',40,'FontWeight','bold','YScale','log');
else
    plot(time(1:nskip:end), ct(1:nskip:end),lineStyle{ls},'LineWidth',5,'Color',colormat(whichline,:));
end
%COLOR SCHEME
%'Color',[0 1 1]);
%'Color',[0 0 1]);


 %   'Color',[0.494117647409439 0.184313729405403 0.556862771511078]);
%'Color',[1 0 1]);

 %   'Color',[0.466666668653488 0.674509823322296 0.18823529779911]);

%'Color',[1 0 0]);
%'Color',[1 1 0]);
%'Color',[0 0 0]);
%'Color',[0 1 0]);


fsave=sprintf('avg%dtraj_%s',ntraj, fname);
%fsave=sprintf('avg20_nbndpair_K%g_set%d.dat',Kd, setNum);
%fsave=sprintf('avg20_nbndpair_box%g_K%g_set%d.dat',box, Kd, setNum);
save(fsave, 'rawDataAvg','-ascii');


%fsave2=sprintf('avg20_meanComplexSize_K%g_dt%g.dat',Kd, deltat);
%save(fsave2, 'mtot','-ascii');
strInterval=sprintf('%g to %g',tstart, tfinish);

[a,b]=find(rawDataAvg(:,2)>1);
maxtime=a(end);
%fit the A(t) data to A(t)=A0Exp(-B0kont)
[BETA1,R,J,COVB,MSE]=nlinfit(rawDataAvg(:,1),rawDataAvg(:,2),@fitcrowd, kon);
CI = nlparci(BETA1,R,'covar',COVB);

%this only fits up until 1 particle is left.
[BETAS1,R,J,COVB,MSE]=nlinfit(rawDataAvg(1:maxtime,1),rawDataAvg(1:maxtime,2),@fitcrowd, kon);
CI = nlparci(BETA1,R,'covar',COVB);

%fit the A(t) data to A(t)=A0Exp(-B0kont)
yfit=log(rawDataAvg(:,2));
t=rawDataAvg(:,1);
[a,b]=find(isinf(yfit));
if(length(a)>0)
    display('CANt take log of zero')
    length(a)
    [a,b]=find(~isinf(yfit));
    t=rawDataAvg(a, 1);
    yfit=log(rawDataAvg(a,2));
end
    
[BETA2,R,J,COVB,MSE]=nlinfit(t, yfit,@fitcrowdlog, kon);
CI2 = nlparci(BETA2,R,'covar',COVB);
[BETAS2,R,J,COVB,MSE]=nlinfit(rawDataAvg(1:maxtime,1),log(rawDataAvg(1:maxtime,2)),@fitcrowdlog, kon);
CI2 = nlparci(BETA2,R,'covar',COVB);

StoreValues=[
time(end),

ntraj,
BETA1/1E6,
CI(1)/1E6,
CI(2)/1E6,
mean(beta)/1e6,
std(beta)/1e6,
std(beta)/sqrt(ntraj)/1e6,
BETA2/1E6,
BETAS1/1e6,
BETAS2/1e6
];
