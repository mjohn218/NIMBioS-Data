%function yhat=fitcrowd(beta, x, N0, V)
function yhat=fitcrowd(beta, x)
N0=100;
V=12500;
yhat=N0*exp(-N0*beta(1)/V*x);