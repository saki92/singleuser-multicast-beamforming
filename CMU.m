function [SNR,SNR_opt] = CMU(R,w,tslots,K,N,sigma)
%CMU
Gk1 = []; Gk2 = [];
epsilon = 0.001; lambda = 1;
for t = 1:tslots
    for k = 1:K
        if w(:,t)'*R(:,:,k)*w(:,t) >= sigma(k,t)
            Gk1 = [Gk1,t];
        else
            Gk2 = [Gk2,t];
        end
        cvx_begin quiet
        variable Rk_est(N,N) complex; %optimization variable
        expression obj; %objective function
        Gk1sum = 0; Gk2sum = 0;
        for n = 1:length(Gk1)
            Gk1sum = Gk1sum + log10(w(:,Gk1(n))'*Rk_est*w(:,Gk1(n))...
                -sigma(k,Gk1(n)));
        end
        for n = 1:length(Gk2)
            Gk2sum = Gk2sum + log10(-(w(:,Gk2(n))'*Rk_est*w(:,Gk2(n)))...
                +sigma(k,Gk2(n)));
        end
        obj = Gk1sum + Gk2sum + logdet(Rk_est);
        maximize obj;
        cvx_end;
        Rcap(:,:,k) = Rk_est;
    end
    wsum = 0;
    for k = 1:K
        wsum = wsum + Rcap(:,:,k)/(w(:,t)'*Rcap(:,:,k)*w(:,t)+epsilon);
    end
    w(:,t+1) = wsum * w(:,t) - (lambda/ceil(0.1*t))*(w(:,1:t)*w(:,1:t)')...
        *w(:,t);
    w(:,t+1) = w(:,t+1)/norm(w(:,t+1));
    for k = 1:K
        sigma(k,t+1) = w(:,t+1)'*Rcap(:,:,k)*w(:,t+1);
        cmuSNR(k) = w(:,t)'*Rcap(:,:,k)*w(:,t);
    end
end
SNR = min(cmuSNR); %minimum SNR among all users
end