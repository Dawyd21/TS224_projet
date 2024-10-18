function [y]=Mon_Welch(x,NFFT,Fe)
% l algorithme du périodogramme de Welch sans chevauchement et sans fenêtre de pondération.
% x est le vecteur contenant les échantillons du signal pour lequel il faut calculer la DSP
% NFFT représente le nombre de points sur lequel les FFT doivent être calculées
% y est l’estimation de la DSP de x
% Fe est la fréquence d’échantillonnage
len=length(x);
k=floor(len/NFFT);
tmp=[];
y=0;
for i=1:k
    tmp=x(1+(i-1)*NFFT:i*NFFT);
    tmp=fftshift(fft(tmp));
    tmp=(abs(tmp).^2)/(NFFT);
    y=y+tmp;
end
y=y/k;
end