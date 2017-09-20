function kt = ktnum(ts,s,D,k)

for i = 1:length(ts)
    t = ts(i);
    integra = integral(@(x)ktnumtobeintegrated(x,t,s,D,k),0,Inf,'RelTol',1e-10,'AbsTol',1e-10);
    kt(i) = k*(1-integra*s*k);
end
