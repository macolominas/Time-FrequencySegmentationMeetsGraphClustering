function [SST1,SST2,omega1,omega2] = my_synchro_order_2(Fg,ft,sigma,Ftg,Ft2g,gamma)
% This function implements Synchrosqueezing up to 2nd order using STFTs
% matrices Ftng.

[neta,nb] = size(Fg);

redun = 1/(ft(2)-ft(1));

omega1 = redun*(repmat(ft',1,nb) + (1/sigma^2)*imag(Ftg./Fg));

N2 = -Fg.*Ftg;
D2 = Ft2g.*Fg - Ftg.^2;
omega2 = redun*(repmat(ft',1,nb) - imag(N2./D2)/(2*pi));


SST1 = zeros(neta,nb);
SST2 = zeros(neta,nb);

for b=1:nb
    for eta=1:neta
        if abs(Fg(eta,b))>gamma
        k = round(omega1(eta,b));
        if k>=1 && k<=neta
            SST1(k,b)  = SST1(k,b) + Fg(eta,b);
        end
        
        k = round(omega2(eta,b));
        if k>=1 && k<=neta
            SST2(k,b)  = SST2(k,b) + Fg(eta,b);
        end

        end
    end
end
        
        
        
        