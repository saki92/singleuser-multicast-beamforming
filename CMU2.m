function [SNR, w, Rcap] = CMU2(R,w,tslots,K,N,sigma)
%CMU
epsilon = 0.0001; lambda = 1;
Gk1 = zeros(K,tslots); Gk2 = zeros(K,tslots);
for t = 1:tslots
    t
    %sigma(:,t)
    for k = 1:K
        if w(:,t)'*R(:,:,k)*w(:,t) >= sigma(k,t)
            Gk1(k,t) = 1;
        else
            Gk2(k,t) = 1;
        end
        cvx_begin quiet
        variable Rk_est(N,N) complex; %optimization variable
        expression obj; %objective function
        Gk1sum = 0; Gk2sum = 0;
        for n = 1:t
            Gk1sum = Gk1sum + Gk1(k,n)*log(real(trace(w(:,n)*w(:,n)'*Rk_est)...
                -sigma(k,n)));
            Gk2sum = Gk2sum + Gk2(k,n)*log(real(-trace(w(:,n)*w(:,n)'*Rk_est)...
                +sigma(k,n))); 
        end
        %obj = Gk1sum + Gk2sum + log_det(0.5*(Rk_est+Rk_est'));
        obj = Gk1sum + Gk2sum + log_det(Rk_est);
        maximize obj;
        cvx_end;
        Rcap(:,:,k) = Rk_est;
        if Rk_est ~= Rk_est
            return;
        end
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
        cmuSNR(k) = w(:,t)'*R(:,:,k)*w(:,t);
    end
    SNR(t) = min(cmuSNR); %minimum SNR among all users
    %w(:,t+1) = w(:,t+1) / SNR(t);
end
end