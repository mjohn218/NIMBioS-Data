function kttot = kt2D(ts,s,D,k)

kt=zeros(length(ts),1);
for i = 1:length(ts)
    t = ts(i);
    integra = integral(@(x)ktnumtobeintegrated(x,t,s,D,k),0,Inf,'RelTol',1e-10,'AbsTol',1e-10);
    kt(i) = k*(1-integra*s*k);
end
kttot=[ts, kt];