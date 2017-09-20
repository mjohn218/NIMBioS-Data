function At = Rev2D_MRE_AB(A0,B0,k,kb,s,D,time)
% Concentrations are in units of molecules per nm^{2}
% Length scale s is in units of nm
% time is in us ({\mu}s)  
% A0(molecules per nm^{2})
% k and D are in units of nm^2/({\mu}s)

scalarshorttime=0.001; %shortime asymptote
scalarlongtime=10; %long time asymptote
%to avoid asymptotis, set shorttime=0 and longtime>time(end)*sig^2/D.

%%% FOR RXN TYPE: A+B<->C
numTsample=1000; %length of the time vector for k(t)
[ktvec,tvec] = ktassister(time,D,k,s,scalarlongtime,scalarshorttime, numTsample);


[~,res] = ode23s(@(t,y) ODEforrxntype4(t,y,k,kb,ktvec,tvec),time,[A0 B0 0]);
At = res(:,1);
Attot=[time, At];

function dy = ODEforrxntype4(t,y,k,kb,ktvec,tvec)

dy = zeros(3,1);    % a column vector
thekt = givekt(t,ktvec,tvec);
dy(1) = -thekt * y(1) * y(2) + kb/k * thekt * y(3);
dy(2) = -thekt * y(1) * y(2) + kb/k * thekt * y(3);
dy(3) =  thekt * y(1) * y(2) - kb/k * thekt * y(3);