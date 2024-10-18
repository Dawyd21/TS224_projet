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
m=2; % parametre pour le lissage de daniell

%% 2.2 - Préambule

variance=3; % sigma
bruit = variance*randn(1,n);

% Fonction d'autocorrelation theorique
R_bruit=cat(2,zeros(1,n-1),ones(1,1));
R_bruit=cat(2,R_bruit,zeros(1,n-1))*variance^2;

% Calcul de la fonction d'autocorrélation
[c_biased, lags] = xcorr(bruit,'biased');
window=[zeros(1,ceil(length(c_biased)/4)),hamming(floor(length(c_biased)/2))',zeros(1,ceil(length(c_biased)/4))];
correlogram=abs(fft(window.*c_biased));
figure,
plot(lags,c_biased);
title("Fonction autocorrelation biaisée");

[c_unbiased, lags] = xcorr(bruit,'unbiased');
figure,plot(freq_axis,c_unbiased)

title("fonction autocorrelation non biaisée");

% Calcul du spectre de puissance
DSP_bruit_WK = fftshift(abs(fft(c_biased))); % avec le theoreme de wiener-khintchine 
DPS_bruit_Welsh= Mon_Welch(bruit,n/100,Fe,n/100-1);% avec periodogramme welsh
DPS_bruit_Barlett= Mon_Barlett(bruit ,n/100,Fe);
DPS_bruit_Daniell= Mon_Daniell(DPS_bruit_Welsh,m);

freq_axis_periodo=linspace(-Fe/2,Fe/2,length(DPS_bruit_Welsh));
figure,hold on
plot(freq_axis_periodo,DPS_bruit_Welsh,LineWidth=2);
plot(freq_axis_periodo,DPS_bruit_Barlett,LineWidth=2);
plot(freq_axis_periodo,DPS_bruit_Daniell,LineWidth=2);
plot(freq_axis_periodo,(variance^2)*ones(1,length(DPS_bruit_Welsh)));
title("DPS du bruit")
legend('Periodogramme de Welsh','Periodogramme de Barlett','Periodogramme de Daniell');

% changement de la taille de l'axe de freq
figure,hold on


plot(freq_axis,DSP_bruit_WK)
plot(freq_axis,correlogram)
plot(freq_axis,ones(1,length(DSP_bruit_WK))*variance^2,LineWidth=2)
title("Comparaison entre la dsp theorique et pratique")
legend('DSP avec WK','correlogramme','theorique')


