%function yhat=fitcrowd(beta, x, N0, V)
function yhat=fitcrowdlog(beta, x)
N0=100;
V=12500;
%lny=lnN-N0/B*beta*t
yhat=log(N0)-N0*beta(1)/V*x;