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
end
for j = 1:K
    cSNR(j) = w'*R(:,:,j)*w;
end
%MU-SLA
%%%Dima%% You need to iterate untill the norm of the difference between the
%%%%% current and the previous beamformer less the 10^-4 or a number of
%%%%% iteration 1000 whatever comes first...As in the paper

tSNR = min(cSNR);
wMU = w;
w = w/sqrt(tSNR);
for it = 1:1
    cvx_begin quiet
    variable w_opt(N, 1) complex; %optimization variable
    expression obj; %objective function
    expression s(K,1) %constraint functions
    for j = 1:K
        s(j) = norm([real(w'*H(:,j)),imag(w'*H(:,j))].')^2+...
            2*[real(w'*H(:,j)),imag(w'*H(:,j))]*([real(w_opt'*H(:,j))...
            ,imag(w_opt'*H(:,j))].'-[real(w'*H(:,j)),imag(w'*H(:,j))].');% Dima%%%%I think the constraint is worng...
        %%please write it as defined in the paper using the square of the norm...It should be a scalar value not a vector
    end
    obj = sum_square_abs(w_opt);
    minimize obj;
    subject to;
    for j = 1:K
        s(j) >= 1;
    end
    cvx_end;
    w = w_opt;
    clear w_opt;
end
w_opt = w;
%Comparing SNRs
w_opt = w_opt/norm(w_opt);  %Normalizing SLA w
%wMU = wMU/norm(wMU); %Normalizing back the MU w again %%%%Dima%%%why you do that you normailzed this to one already!!!
for j = 1:K
    muSNR(j) = (wMU'*R(:,:,j)*wMU)/sigma; %SNR with MU w
    muSNR_opt(j) = (w_opt'*R(:,:,j)*w_opt)/sigma; %SNR with MU-SLA w
end
SNR = min(muSNR); SNR_opt = min(muSNR_opt); %minimum SNR among all users
%SNR = muSNR; SNR_opt = muSNR_opt;
end