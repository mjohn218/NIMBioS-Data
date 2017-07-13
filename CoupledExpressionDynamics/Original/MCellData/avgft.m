function [avgA, avgR] = avgft(numseeds, filenameA, filenameR)
%Usage: plotft(24, 'A.Cube.dat', 'R.Cube.dat')

numfol = numseeds;
avglengthA = 0;
avglengthR = 0;

%Extract Time
time = dlmread('seed_00001/A.Cube.dat');
time = time(:,1);

for i = 1:numfol
    s = sprintf('%2.2d',i);
    A = dlmread(['seed_000' s  '/' filenameA]);
    A = A(:,2);
    R = dlmread(['seed_000' s '/' filenameR]);
    R = R(:,2);
    
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

    difft = diff(time);
    Fs = 1/difft(1);
    f = Fs*(0:(L/2))/L;

    [y,i]=max(P1(2:end));

    avglengthA = avglengthA + 1/f(i+1);
    
    % Find wavelength of R
    L = length(R);
    if mod(L,2) ~= 0
        R = R(1:end-1,:);
    end

    Y = fft(R);
    L = length(R);

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    difft = diff(time);
    Fs = 1/difft(1);
    f = Fs*(0:(L/2))/L;

    [y,i]=max(P1(2:end));

    avglengthR = avglengthR + 1/f(i+1);
    
end
avgA = avglengthA/numseeds;
avgR = avglengthR/numseeds; 
