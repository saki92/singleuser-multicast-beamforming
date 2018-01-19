clear;
K = [25,100,200,300,400,500]; %no of receivers
N = [20,40,60,80,100]; %no of tx antennas
iteration = 200; %no of iterations for MU
%w = randn(N,1) + 1i*randn(N,1); %initial w vector
sigma = 0.1;
RicSig = 0.1;
G = 25; %no of groups for rician channel
Rician = 1;
%for k = 1:length(K) %for all users
for k = 1:length(N) %for all antennas
    %Np = 20;
    Np = N(k)
    %Kp = K(k)
    Kp = 450;
    w = randn(Np,1) + 1i*randn(Np,1);
    for m = 1:100 %for averaging over 100 channel instances
        if Rician == 1
            Ht = 1/sqrt(2)*(randn(Np,G)+1i*randn(Np,G));
            Hg = [repmat(Ht,1,fix(Kp/G)),Ht(:,1:mod(Kp,G))];
            Hadd = RicSig*(1/sqrt(2)*(randn(Np,Kp)+1i*randn(Np,Kp)));
            H = Hg + Hadd;
        else
            H = 1/sqrt(2)*(randn(Np,Kp)+1i*randn(Np,Kp));
        end
        [SNR(m,:), SNR_opt(m,:)] = MUSLA(H,w,iteration,Kp,Np,sigma);
    end
    amSNR(k,:) = mean(SNR,1); amSNR_opt(k,:) = mean(SNR_opt,1);
end
%plot(K,10*log10(real(amSNR_opt)),'b-o'); hold on
%plot(K,10*log10(real(amSNR)),'-xk');
plot(N,10*log10(real(amSNR_opt)),'b-o'); hold on
plot(N,10*log10(real(amSNR)),'-xk');
%plot(K,20*log10(real(amSNR_opt)),K,20*log10(real(amSNR)));
%xlabel('Number of Users');
xlabel('Number of Antennas');
ylabel('Average minimum SNR(dB)');
legend('MU-SLA','MU')