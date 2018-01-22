function SNR = CMU(R,w,tslots,K,N,sigma)
%CMU
epsilon = 0.01; lambda = 1;
Gk1 = cell(K,1); Gk2 = cell(K,1);
for t = 1:tslots
    t
    %sigma(:,t)
    for k = 1:K
        if w(:,t)'*R(:,:,k)*w(:,t) >= sigma(k,t)
            Gk1{k} = [Gk1{k},t];
        else
            Gk2{k} = [Gk2{k},t];
        end
        cvx_begin
        variable Rk_est(N,N) complex; %optimization variable
        expression obj; %objective function
        Gk1sum = 0; Gk2sum = 0;
        for n = 1:length(Gk1{k})
            Gk1sum = Gk1sum + log(real(trace(w(:,Gk1{k}(n))*w(:,Gk1{k}(n))'*Rk_est)...
                -sigma(k,Gk1{k}(n))));
        end
        for n = 1:length(Gk2{k})
            Gk2sum = Gk2sum + log(real(-trace(w(:,Gk2{k}(n))*w(:,Gk2{k}(n))'*Rk_est)...
                +sigma(k,Gk2{k}(n))));
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
        cmuSNR(k) = w(:,t)'*Rcap(:,:,k)*w(:,t);
    end
    SNR(t) = min(cmuSNR); %minimum SNR among all users
end
end