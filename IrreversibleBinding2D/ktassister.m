function [ktvec,tvec] = ktassister(time,D,k,s,scalarlongtime,scalarshorttime, numsample)
%creates sampled k(t)

%numsample = 1000;

ttemp = logspace(log10(time(2)),log10(time(end)),numsample);
ttemp(1) = time(2);
ttemp(end) = time(end);
tvec = [0;ttemp(:)];
ktvec = createkt(tvec,D,k,s,scalarlongtime,scalarshorttime);