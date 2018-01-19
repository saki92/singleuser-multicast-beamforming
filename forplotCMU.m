clear; clc;
K = 10;
N = 2;
tslots = 400;
w = randn(N,1) + 1i*randn(N,1); %initial w vector
for m = 1:100 %for averaging over 100 channel instances
    m
    sigma = randn(K,1);
    for k = 1:K
        M = 1/sqrt(2)*(randn(N,N)+1i*randn(N,N));
        R(:,:,k) = M*M';
        sigma(k) = w'*R(:,:,k)*w;
    end
    SNR(m,:) = CMU(R,w,tslots,K,N,sigma);
end
amSNR = mean(SNR,1);
plot(1:tslots,10*log10(real(amSNR)),'b-o');
xlabel('Number of time slots');
ylabel('Average minimum SNR(dB)');
%legend('MU-SLA','MU')