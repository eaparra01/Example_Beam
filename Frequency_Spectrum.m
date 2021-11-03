function [f,P1] = Frequency_Spectrum(time,data)

    Fs = 1/(time(364+1)-time(364));
    T = 1/Fs;

    L = length(data);
    Y = fft(data);

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    f = Fs*(0:(L/2))/L;
