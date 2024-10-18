function [y]=Mon_Welch(x,NFFT,Fe,overlap)
% l algorithme du périodogramme de Welch sans chevauchement et sans fenêtre de pondération.
% x est le vecteur contenant les échantillons du signal pour lequel il faut calculer la DSP
% NFFT représente le nombre de points sur lequel les FFT doivent être calculées
% y est l’estimation de la DSP de x
% Fe est la fréquence d’échantillonnage
len=length(x);
ptr = NFFT - overlap;   % decalage
k=floor((len-overlap)/ptr); % nbr d'intervalle
tmp=[];
y=0;
window=hamming(NFFT);
U = sum(window.^2) / NFFT; % facteur de correction
for i=1:k
    if(i*ptr+NFFT<=len)
        tmp=x(1+(i-1)*ptr:(i-1)*ptr+NFFT);
        tmp=tmp.*transpose(window);
        tmp=fftshift(fft(tmp));
        tmp=(abs(tmp).^2)/(NFFT*Fe*U);
        y=y+tmp;
    end
end
y=y/k;
end