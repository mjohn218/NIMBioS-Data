function plotft(dattype, time,dat)
%Usage: plotft('A', timeseriesA,timeseries dataA)
%Usage: plotft('R', timeseriesR,timeseries dataR)

L = length(dat);
if mod(L,2) ~= 0
    dat = dat(1:end-1,:);
end

Y = fft(dat);
L = length(dat);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

difft = diff(time);
Fs = 1/difft(1);
f = Fs*(0:(L/2))/L;

plot(f,P1);
title(['Single-Sided Amplitude Spectrum of ' dattype '(t)'])
xlabel('f (Hz)')
ylabel(['|' dattype '(f)|'])

[y,i]=max(P1(2:end));

disp(['Predicted wavelength: ' num2str(1/f(i+1)) ' s.'])