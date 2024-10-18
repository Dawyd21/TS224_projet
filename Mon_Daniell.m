function [x] = Mon_Daniell(x, m)
% x est le periodogramme a lissé
% m est le parametre tel qu'on considere les valeurs entre t-m et t+m
% inclus, on note nbr le nombre de valeurs considere
len=length(x);
nbr=2*m+1;
for i=1:len
    if i-m>1 && i+m<len
        tmp=x(i);
        for k=1:m
            tmp=tmp+x(i-k)+x(i+k);
        end
        x(i)=tmp/nbr;
    end  
end

end