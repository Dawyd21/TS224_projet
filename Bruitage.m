function [x] = Bruitage(x,SNR)
% x est le signal a bruite
% on rajoute du bruit par dessus
% on choisi d'ajouter une BBGC avec randn
% RSB=10*log(P_s/P_b)
% Puissance_bruit=Puisance_signal/(10^(SNR/10))
% E_b= E_s/(10^(SNR/10))
len=length(x);
E_s=mean(x.*x);
sigma2=E_s/(10^(SNR/10));
bruit=sqrt(sigma2)*randn(1,len);
x=x+bruit;
end
