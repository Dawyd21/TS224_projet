%% Projet TS224
% Ce script génère un bruit blanc gaussien, calcule et compare ses
% autocorrélations (théorique, biaisée, et non biaisée) et analyse son
% spectre de puissance à l'aide de la transformée de Fourier et de la
% méthode de Welch.

%% Initialisation
clc
clear
close all

%% Paramètres
n = 10000;  % Nombre de points du signal
Fe = 8000;  % Fréquence d'échantillonnage
time_axis = 0:1/Fe:(1/Fe)*n - 1/Fe;  % Axe temporel
freq_axis = -Fe * n + Fe : Fe : Fe * n - Fe;  % Axe fréquentiel

%% 2.2 - Préambule

% Paramètres du bruit
variance = 3; % Variance (sigma^2)
bruit = variance * randn(1, n);  % Génération du bruit blanc gaussien centré

% Fonction d'autocorrélation théorique
% Pour un bruit blanc, l'autocorrélation est nulle sauf à l'origine où elle
% vaut la variance (ici, 3^2 = 9).
R_bruit = cat(2, zeros(1, n - 1), ones(1, 1));
R_bruit = cat(2, R_bruit, zeros(1, n - 1)) * variance^2;

% Affichage de l'autocorrélation théorique
figure;
stem(-n + 1:n - 1, R_bruit); % Utilisation de stem pour mieux visualiser le pic
title('Autocorrélation théorique du bruit blanc');
xlabel('Décalage');
ylabel('Autocorrélation');
grid on;

% Calcul de la fonction d'autocorrélation estimée (biaisée et non biaisée)
[c_biased, lags] = xcorr(bruit, 'biased');
[c_unbiased, lags] = xcorr(bruit, 'unbiased');

% Affichage des autocorrélations estimées
figure;
plot(lags, c_biased);
title('Fonction autocorrélation biaisée');
xlabel('Décalage');
ylabel('Autocorrélation');
grid on;

figure;
plot(lags, c_unbiased);
title('Fonction autocorrélation non biaisée');
xlabel('Décalage');
ylabel('Autocorrélation');
grid on;

% Analyse des résultats :
% On observe que l'autocorrélation théorique a un pic à l'origine (décalage
% zéro) égal à la variance et est nulle ailleurs. L'autocorrélation biaisée
% et non biaisée diffèrent par la manière dont elles traitent les valeurs
% proches des bords (le biais augmente lorsque le décalage augmente).

% Calcul du spectre de puissance à l'aide de la transformée de Fourier
% (avec le théorème de Wiener-Khinchine)
DSP_bruit = fftshift(abs(fft(c_biased, 2 * n))); 

% Calcul de la DSP avec la méthode de Welch
DPS_bruit = Mon_Welch(bruit, n / 100, Fe); % Assurez-vous que la fonction Mon_Welch est correcte

% Affichage des spectres de puissance
figure;
semilogy(DPS_bruit);
title('DPS du bruit (méthode de Welch)');
xlabel('Fréquence (Hz)');
ylabel('DPS');
grid on;

% Comparaison entre la DSP théorique (constante) et celle obtenue
freq_axis = -Fe * n + Fe : Fe : Fe * n;
figure;
hold on;
plot(freq_axis, DSP_bruit);
plot(freq_axis, ones(1, length(DSP_bruit)) * variance^2, 'r', 'LineWidth', 2);
title('Comparaison entre la DSP théorique et pratique');
xlabel('Fréquence (Hz)');
ylabel('DPS');
legend('Pratique', 'Théorique');
grid on;

% Analyse des résultats :
% La densité spectrale de puissance théorique pour un bruit blanc est
% constante et égale à la variance (ici, 9). La DSP obtenue à partir de la
% transformée de Fourier et celle obtenue par la méthode de Welch devraient
% approcher cette valeur. Les écarts peuvent être dus à la nature aléatoire
% du bruit et à la taille finie de l'échantillon.
