function [y, f] = Mon_Daniell(x, Nfft, Fe)
    X = fft(x, Nfft);
    Px = abs(X).^2 / Nfft;
    L = 5; 
    window = ones(1, L)/L; 

    Px_lisse = conv(Px, window, 'same'); 

    f = (0:(Nfft-1)) * (Fe/Nfft);
    f = f - Fe/2;

    y = fftshift(Px_lisse);
end