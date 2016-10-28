close all
clear;
clc;

max = 5e4;
cnt = zeros(max,13);
divisor = cnt;
maxfol = 24;
% 
% V = 0.02;
% NA = 6.022e23;
% umtodm = 1e-5;
% Vindm3 = V*umtodm^3;
% 
% factor = 1/Vindm3/NA;
cd('gonder')
for i = 1:maxfol
   
    cd(num2str(i))
    mat = dlmread('timestat.dat');    
    mat = mat(1:2:end,:);
    
    mat = mat(1:max,:);
    
    [nr,nc]=size(mat);
    if nc<13
        mat = [mat zeros(nr,13-nc)];
    end
    dt = mat(2,1)-mat(1,1);
    t = 0:dt:(max-1)*dt;
    
%     divisor = divisor + [ones(max,1) [ones(nr,12); zeros(max-nr,12)]];
    
%     mat = [t(:) [mat(:,2:end); zeros(max-nr,12)]];

    %     mat = mat(1:max,:);
    cnt = cnt + mat;
    cd('..')
%     eqs(i) = mat(end,7);
end
cd('..')

mat = cnt/maxfol;

dlmwrite('A.dat',[mat(:,1) mat(:,3)],'delimiter','\t','precision',6);

% disp([mean(eqs) std(eqs)/sqrt(maxfol)])
% 
% M(1,:) = [mean(eqs) std(eqs)/sqrt(maxfol)];
% 
% cnt = cnt./divisor; %maxfol;
% cnt(1,:) = [0,35631,120,120,0,0,0,0,0,0,0,0,0];
% % pmat = [cnt(:,3) cnt(:,4) cnt(:,6) cnt(:,7) cnt(:,5) cnt(:,8)]/1927;
% 
% % pmat = [cnt(:,3) cnt(:,4) cnt(:,6) cnt(:,7) cnt(1,3)-(cnt(:,3)+cnt(:,6)+cnt(:,12)) cnt(:,9)+cnt(:,10)+cnt(:,13)];
% 
% sum = (cnt(:,5)+cnt(:,6)+cnt(:,7))/2;
% 
% semilogx(cnt(:,1),[cnt(:,3) cnt(:,4) sum-cnt(:,7) sum-cnt(:,6) sum-cnt(:,5) cnt(:,13)+cnt(:,10)]/L/L,'LineWidth',2);
% legend('Free A Cyt.','Free B Cyt.','Free A Mem.','Free B Mem.','AB Cyt.','AB Mem.','Location','SouthWest');
% % title(['k_a:0.1 nm^3/{\mu}s AP2:' num2str(round(cnt(1,3))) ' EPN:' num2str(round(cnt(1,4))) ' PIP:' num2str(round(cnt(1,2))) ' {\sigma}: 5nm'],'FontSize',16)
% xlabel('Time (s)','FontSize',16);
% set(gca,'FontSize',16);
% % ylabel('Concentration (nm^{-2})','FontSize',16);
% ylabel('#Molecules','FontSize',16);
% set(gcf,'color','w');
% % 
% % hold on
% % xx = xlim;
% % mat = dlmread('plot.dat');
% % semilogx(mat(:,1),mat(:,2:end)/1927,'k--')
% % % ylim([0 5000])
% % xlim([xx(1) xx(2)])
% % saveas(gcf, 'ka0.02.fig', 'fig')
% 
% ylim([0 140]);
% 
% xlim([5e-5 1e-1])
% 
% te = input('Last time');
% %                1        2         3             4            5               6
% lastdata = [cnt(te,3) cnt(te,4) sum(te)-cnt(te,7) sum(te)-cnt(te,6) sum(te)-cnt(te,5) cnt(te,13)+cnt(te,10)];
% 
% sol1 = lastdata(1)+lastdata(3);
% sol2 = lastdata(2)+lastdata(4);
% comp = lastdata(5)+lastdata(6);
% 
% Kae = comp/sol1/(sol2*factor);
% 
% disp(Kae)
% 
% % x=cnt(1:te,1);
% % v=cnt(1:te,13)+cnt(1:te,10);
% % 
% % vextrap = interp1(x,v,cnt(1:10:end,1),'pchip','extrap');
% % hold on
% % semilogx(cnt(1:10:end,1),vextrap,'k--');
% 
% 
% % index = round(logspace(0,floor(log10(max)),40));
% % index = unique([index(:);length(cnt)]);
% % 
% % complexes = cnt(index,7);
% % time = cnt(index,1);
% % 
% % % complexes = cnt(:,7);
% % % time = cnt(:,1);
% % Rct = (complexes-complexes(end))/(complexes(1)-complexes(end));
% % Rctdiff = diff(Rct);
% % Rctdiff = Rctdiff(:);
% % timeedited = 0.5 * (time(1:end-1)+time(2:end));
% % timeedited = timeedited(:);
% % equiltime = -trapz(timeedited, timeedited.*Rctdiff);
% % eqtimesq = -trapz(timeedited, timeedited.*timeedited.*Rctdiff);
% % eqvar = eqtimesq - equiltime*equiltime;
% % eqvarnodim = eqvar/(equiltime*equiltime);
% 
% index = round(logspace(0,floor(log10(max)),40));
% index = unique([index(:);length(cnt)]);
% 
% complexes = cnt(index,7);
% time = cnt(index,1);
% 
% complexes = pchip(time,complexes,cnt(:,1));
% time = cnt(:,1);
% 
% Rct = (complexes-complexes(end))/(complexes(1)-complexes(end));
% F = 1-Rct;
% 
% equiltime = time(end)*F(end)-time(1)*F(1)-trapz(time,F);
% eqtimesq = F(end)*time(end)^2-F(1)*time(1)^2-2*(time(end)*F(end)-time(1)*F(1))+2*trapz(time,F);
% eqvar = eqtimesq - equiltime*equiltime;
% 
% disp([equiltime eqvar])
% M(2,:) = [equiltime eqvar];
% dlmwrite('file.txt',M,'delimiter','\t','precision','%.6f')
% % hold on
% % scatter(time,complexes/L/L,'*k')
% % set(gca, 'xscale', 'log')