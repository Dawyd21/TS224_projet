%% Projet TS224
% Ce script génère un bruit blanc gaussien, calcule et compare ses
% autocorrélations (théorique, biaisée, et non biaisée) et analyse son
% spectre de puissance à l'aide de la transformée de Fourier et de la
% méthode de Welch.

%% Initialisation
clc
clear
close all

%% Parametre
load fcno03fz.mat
load data_Weierstrass.mat
n=length(fcno03fz);
Fe=8000;
time_axis=0:1/Fe:(1/Fe)*n-1/Fe;
freq_axis=-Fe*n+Fe:Fe:Fe*n-Fe;
m=2; % parametre pour le lissage de daniell

% RSB definition a=sqrt(P_signal/P_bruit*10**(-RSB/10))


%% 2.2 - Préambule

% Paramètres du bruit
variance = 3; % Variance (sigma^2)
bruit = variance * randn(1, n);  % Génération du bruit blanc gaussien centré

% Fonction d'autocorrélation théorique
% Pour un bruit blanc, l'autocorrélation est nulle sauf à l'origine où elle
% vaut la variance (ici, 3^2 = 9).
R_bruit = cat(2, zeros(1, n - 1), ones(1, 1));
R_bruit = cat(2, R_bruit, zeros(1, n - 1)) * variance^2;


% Calcul de la fonction d'autocorrélation
[c_biased, lags] = xcorr(bruit,'biased');
window=[zeros(1,ceil(length(c_biased)/4)),hamming(floor(length(c_biased)/2))',zeros(1,ceil(length(c_biased)/4))];
correlogram=abs(fft(window.*c_biased));
figure,
plot(lags,c_biased);
title("Fonction autocorrelation biaisée");

[c_unbiased, lags] = xcorr(bruit,'unbiased');
figure,plot(freq_axis,c_unbiased)

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

%% bruitage du signal de parole

signal_RSB10=Bruitage(fcno03fz',10);
signal_weier=Bruitage(data{1,1}',10);

%% Spectre du signal de parole
[~,F,T,P_bruit]=spectrogram(signal_RSB10,256,128,256);
[~,~,~,P_original]=spectrogram(fcno03fz,256,128,256);
P_bruit_log=10*log10(P_bruit);
P_original_log=10*log10(P_original);



figure,
subplot(211)
imagesc(T, F, P_bruit_log)
axis xy 
xlabel('Temps (s)')
ylabel('Fréquence Normalisée')
colorbar
colormap jet
subplot(212)
plot(time_axis,signal_RSB10)
title("Signal de parole avec un bruit RSB=10")
figure,
subplot(211)
imagesc(T, F, P_original_log)
axis xy
xlabel('Temps (s)')
ylabel('Fréquence (Hz)')
colorbar
colormap jet
subplot(212)
plot(time_axis,fcno03fz)
title("Signal de parole original")


%% DFA

% --------------------------creation du signal test----------------------
variance_test=3; % sigma
esperance_test=10;
M=1000;
L=5;
N=floor(M/L);
signal_test = variance_test*randn(1,M)+esperance_test;
axe=0:M-1;

% --------------------------standardisation du signal--------------------
signal_test_centre=(signal_test-mean(signal_test))/std(signal_test);
%integral_signal_centre=trapz(signal_test_centre);

% --------------------------decoupe du signal en N morceaux--------------
signal_decoupe=[];
for i=1:N
    signal_decoupe=[signal_decoupe ,signal_test_centre((i-1)*L+1:i*L)'];
end
% --------------------------OLS de chaque morceau------------------------
Beta=zeros(2,N);
for i=1:N
    X=[ones(L,1),axe((i-1)*L+1:i*L)'];
    X_t=X';
    Beta(:,i)=inv(X_t*X)*X_t*signal_decoupe(:,i);
end
figure,plot(axe,signal_test_centre,LineWidth=2)
hold on
signal_OLS=zeros(1,M);
for i=1:N
    signal_OLS(1,(i-1)*L+1:i*L)=Beta(1,i)+axe((i-1)*L+1:i*L)*Beta(2,i);
end
plot(axe,signal_OLS,LineWidth=2)

% residu
residu=signal_test_centre-signal_OLS;
plot(axe,residu,LineWidth=2)
legend("Signal centré","signal OLS","residu")
hold off

% carre de la fonction de fluctuation
F_2_test=mean(residu.*residu);
%% DFA d'un signal
signal=fcno03fz;
len=floor(length(fcno03fz)/2);
axe=1:len;
F_2=zeros(1,len);
for i=2:len
    F_2(i-1)=DFA(fcno03fz',i,1);
end
% tracer de F_2 et estimation de alpha
alpha=zeros(2,1);
X_log=[ones(len,1),log(axe(1:len))'];
alpha(:)=X_log\F_2';
figure,semilogy(log(axe),log(F_2))
title("Log(F_2) en fonction de log(x)")
disp("H est egale a")
disp(alpha(2,1)-1)

%% poly 2nd degre
len=floor(length(fcno03fz)/3);
axe_poly2=1:len;
F_2_poly2=zeros(1,len);
for i=3:len
    F_2_poly2(i-2)=DFA(fcno03fz',i,2);
end

% tracer de F_2 et estimation de alpha
alpha_poly2=zeros(2,1);
X_log_poly2=[ones(len,1),log(axe_poly2(1:len))'];
alpha_poly2(:)=X_log_poly2\F_2_poly2';
figure,semilogy(log(axe_poly2),log(F_2_poly2))
title("Log(F_2) en fonction de log(x) avec poly 2nd degre")
disp("H est egale a")
disp(alpha(2,1)-1)


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
