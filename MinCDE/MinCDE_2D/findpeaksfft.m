function wl = findpeaksfft(A,distance)

r = 10;
y = interp(distance,r); % increase sampling frequency by r
y(y>max(distance)) = [];
y(1) = 0;

A = interp1(distance,A,y);
A = A-mean(A);
distance = y;

% Find wavelength of A
    L = length(A);
    if mod(L,2) ~= 0
        A = A(1:end-1,:);
    end

    Y = fft(A);
    L = length(A);

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    difft = diff(distance);
    Fs = 1/difft(1);
    f = Fs*(0:(L/2))/L;

    [~,i]=max(P1);

    wl = 1/f(i);
    
    subplot(1,2,1)
    plot(distance(1:length(A)),A)
    xlabel('Length (um)', 'fontsize',14);
    ylabel('EminDT');
    subplot(1,2,2)
    plot(1./f,P1)
    xlabel('Wavelength Histogram (um)', 'fontsize',14);
    ylabel('EminDT');
    title(['Wavelength ' num2str(wl)])
    pause(1)
    