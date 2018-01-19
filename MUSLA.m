function [SNR,SNR_opt] = MUSLA(H,w,iteration,K,N,sigma)
%MU
%%%Dima%% You need to iterate untill the norm of the difference between the
%%%%% current and the previous beamformer less the 10^-4 or a number of
%%%%% iteration 1000 whatever comes first...As in the paper
for j = 1:K
    R(:,:,j) = H(:,j)*H(:,j)';  %Channel cor-matrix for all users
end
for i = 1:iteration
    summ = zeros(N,N);
    for j = 1:K
        summ = summ + (R(:,:,j)/(w'*R(:,:,j)*w + 10e-9));
    end
    w_til = summ * w;
    w = w_til / norm(w_til);
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