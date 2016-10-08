load('time_osci.mat');
load('ode_osci_a.mat');
load('ode_osci_a.mat');
load('smol_osci_a.mat');
load('smol_osci_r.mat');
load('mcell_osci_a.mat');
load('mcell_osci_r.mat');

figure(1);
subplot(2,1,1);
hold on;
odeplot = plot(excelt_osci, ode_osci_a, 'k');
smolplot = plot(excelt_osci, smol_oci_a, 'g');
mcellplot = plot(excelt_osci, mcelldat_a, 'b');
title('Number of A Molecules');
xlabel('time(s)');
ylabel('N_{A}(t)');
legend([odeplot,smolplot,mcellplot],'ODE','Smoldyn','MCell');

subplot(2,1,2);
hold on;
plot(excelt_osci, ode_osci_r, 'k');
plot(excelt_osci, smol_osci_r, 'g');
plot(excelt_osci, mcelldat_r, 'b');
title('Number of R Molecules');
xlabel('time(s)');
ylabel('N_{R}(t)');
