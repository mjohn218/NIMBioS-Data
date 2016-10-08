function rigid()

options = odeset('RelTol',1e-4,'AbsTol',1e-5);
[T,Y] = ode45(@rigid2,[0 0.02],0,options);

A = 4.835975716; %um2
A0 = 602;
AofT = A0 - Y*A;

dlmwrite('A.dat',[T AofT],'delimiter','\t','precision',6);

semilogx(T*1e6,AofT);
hold on
load FPR.dat
semilogx(FPR(:,1)*1e6,FPR(:,2),'LineWidth',3)

function dy = rigid2(t,y)

Kf = 84.21; % s-1*uM-1
dy = 0;    % a column vector
NA = 6.022e23;
B0 = 6045; %molecules/um2
A0 = 1; %uM

V = 1; %1um3
A = 4.835975716; %um2
Vindm3 = V*(1e-5)^3;

Aleftover = A0-1e6*y/NA*A/Vindm3; %in uM
Bleftover = B0-y; %in molecules/um2

dy(1) = Kf * Bleftover * Aleftover; %molecules/um2
