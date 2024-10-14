%% Projet TS224
%% Initialisation
clc
clear
close all
%% Constantes
n=10000;
Fe=8000;
time_axis=0:1/Fe:(1/Fe)*n-1/Fe;
freq_axis=-Fe*n+Fe:Fe:Fe*n-Fe;

%% 2.2

variance=3; % sigma
bruit = variance*randn(1,n);

% Fonction d'autocorrelation theorique
R_bruit=cat(2,zeros(1,n-1),ones(1,1));
R_bruit=cat(2,R_bruit,zeros(1,n-1))*variance^2;

% Calcul de la fonction d'autocorrélation
[c_biased, lags] = xcorr(bruit,'biased');
[c_unbiased, lags] = xcorr(bruit,'unbiased');

figure,plot(freq_axis,c_biased);
title("Fonction autocorrelation biaisée");

figure,plot(freq_axis,c_unbiased)

title("fonction autocorrelation non biaisée");

% Calcul du spectre de puissance
DSP_bruit = fftshift(abs(fft(c_biased, 2*n))); % avec le theoreme de wiener-khintchine 
DPS_bruit=Mon_Welch(bruit,n/100,Fe);% avec periodogramme welsh

figure,semilogy(DPS_bruit)
title("DPS du bruit")

freq_axis=-Fe*n+Fe:Fe:Fe*n;
figure,hold on
plot(freq_axis,DSP_bruit)
plot(freq_axis,ones(1,length(DSP_bruit))*variance^2,Color="red")
title("Comparaison entre la dsp theorique et pratique")
legend('pratique','theorique')


