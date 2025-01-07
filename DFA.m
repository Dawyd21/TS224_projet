function F_2=DFA(signal,L,puissance)

% Parametres
M=length(signal);
N=floor(M/L);
axe=0:M-1;

% ----standardisation du signal----------------
signal_centre=(signal-mean(signal))/std(signal);

% ---------decoupe du signal en N morceaux--------------
signal_decoupe=zeros(L,N);
for i=1:N
    signal_decoupe(:,i)=signal_centre((i-1)*L+1:i*L)';
end

% ---------OLS de chaque morceau------------------------

if(puissance==1)
    Beta=zeros(2,N);
    for i=1:N
        X=[ones(L,1),axe((i-1)*L+1:i*L)'];
        Beta(:,i)=X\signal_decoupe(:,i);
    end
elseif(puissance==2)
    Beta=zeros(3,N);
    for i=1:N
        %X=[ones(L,1),axe((i-1)*L+1:i*L)',axe((i-1)*L+1:i*L)'.^2];
        %X_t=X';
        %Beta(:,i)=(X_t*X)\(X_t*signal_decoupe(:,i));
        X=[ones(L,1),axe((i-1)*L+1:i*L)',axe((i-1)*L+1:i*L)'.^2];
        Beta(:,i)=X\signal_decoupe(:,i);
    end
end

signal_OLS=zeros(1,M);

if(puissance==1)
    for i=1:N
        signal_OLS(1,(i-1)*L+1:i*L)=Beta(1,i)+axe((i-1)*L+1:i*L)*Beta(2,i);
    end
elseif(puissance==2)
    for i=1:N
        signal_OLS(1,(i-1)*L+1:i*L)=Beta(1,i)+axe((i-1)*L+1:i*L)*Beta(2,i)+axe((i-1)*L+1:i*L).^2*Beta(3,i);
    end
end

% residu
residu=signal_centre-signal_OLS;

% carre de la fonction de fluctuation

F_2=mean(residu.*residu);

end