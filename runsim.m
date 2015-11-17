%% Initialize Parameters

clc
clear

alpha = 0;    %coefficient for block fading model
beta = 0.8^2;  % Attenuation loss from non-direct antennas
n0 = 10^(-3);    %noise variance

Nt = 2;  %Nt antennas for each transmitter
Nr = 2;  %Nr antennas for each receiver
M = 3;   %number of users
Pmax = ones(1, M);    %maximum power for each user
upower = ones(1, M);   %power for unicastmpower 
mpower = ones(1, M);   %power for multicast
upower = sqrt(upower); % Change power to voltage
mpower = sqrt(mpower); % Change power to voltage

iternums = 1:10; % number of iterations
N_realization = 500; % Number of times to run simulation
traininglength = 20; % traininglength 2M

averagerateu = zeros(N_realization, length(iternums));
averageratem = zeros(N_realization, length(iternums));
averagerateu_w = zeros(N_realization, length(iternums));
averageratem_w = zeros(N_realization, length(iternums));

%% Start Loop
for realization_idx = 1 : N_realization
        realization_idx
    H = zeros(Nr,Nt,M,M); % forward channel
    Z = zeros(Nr,Nt,M,M); % reciprocal(backward) channel
    for Tr_user_idx = 1:M
        for Re_user_idx = 1:M
            if Re_user_idx==Tr_user_idx
                H(:,:,Re_user_idx,Tr_user_idx) = (randn(Nr,Nt)+1i*randn(Nr,Nt))/sqrt(2);  %~CN(0,1) %Create channels for each user
                Z(:,:,Tr_user_idx,Re_user_idx) = H(:,:,Re_user_idx,Tr_user_idx)';
            else
                H(:,:,Re_user_idx,Tr_user_idx) = (randn(Nr,Nt)+1i*randn(Nr,Nt))/sqrt(2/beta); %Create interference channels
                Z(:,:,Tr_user_idx,Re_user_idx) = H(:,:,Re_user_idx,Tr_user_idx)';
            end
        end
    end  
    
    M1 = traininglength/2;
    M2 = traininglength/2;
    
    %% one iteration per block
   
    InitialGu = rand(Nr, M) + 1i*rand(Nr, M);    % receive filter unicast
    InitialGm = rand(Nr, M) + 1i*rand(Nr, M);    % receive filter multicast
    for k = 1:M
        InitialGu(:, k) = InitialGu(:, k)./norm(InitialGu(:, k));
        InitialGm(:, k) = InitialGm(:, k)./norm(InitialGm(:, k));
    end
    
    Gu = InitialGu;
    Gm = InitialGm;
    Gu_w = InitialGu;
    Gm_w = InitialGm;
    Vu = zeros(Nt, M); % beamformer unicast
    Vm = zeros(Nt, M); % beamformer multicast
    Vu_w = zeros(Nt, M); % beamformer unicast
    Vm_w = zeros(Nt, M); % beamformer multicast
    
    for numiters = 1:length(iternums)
        BfwBr = sign(randn(M1,M));    %training symbols at the transmitter broadcast
        BbwBr = sign(randn(M1,M));    %training symbols at the transmitter broadcast
        for i = 2:M
        BfwBr(:, i) = BfwBr(:, 1);
        BbwBr(:, i) = BbwBr(:, 1);
        end
        Bfw = sign(randn(M1,M));    %training symbols at the transmitter unicast
        Bbw = sign(randn(M2,M));   %training symbols at the receiver unicast
        
        %% bi-directional training
            %%LS algorithm
            %%phase 1: backward training to update beamformer
            [Vu, Vm] = LS_backward(Z, Gu, Gm, M2, n0, Bbw, BbwBr, upower, mpower);
            [Vu_w, Vm_w] = Wiener_backward(Z, Gu, Gm, M2, n0, Bbw, BbwBr, upower, mpower);
            
            %%phase 2: forward training to update receive filter
            [Gu, Gm] = LS_forward(H, Vu, Vm, M1, n0, Bfw, BfwBr, upower, mpower);
            [Gu_w, Gm_w] = Wiener_forward(H, Vu, Vm, M1, n0, Bfw, BfwBr, upower, mpower);
            
        averagerateu(realization_idx, numiters) = calculate_rateu(H, n0, Vu, Gu, Vm, upower, mpower);
        averageratem(realization_idx, numiters) = calculate_ratem(H, n0, Vm, Gm, Vu, upower, mpower);
        averagerateu_w(realization_idx, numiters) = calculate_rateu(H, n0, Vu_w, Gu_w, Vm_w, upower, mpower);
        averageratem_w(realization_idx, numiters) = calculate_ratem(H, n0, Vm_w, Gm_w, Vu_w, upower, mpower);
    end
            
    
end

hold on
plot(iternums, mean(averagerateu)+mean(averageratem), 'b',iternums, mean(averagerateu_w)+mean(averageratem_w), 'k');
xlabel('Number of iterations')
ylabel('Rates')
title('Rates vs number of iterations at 0dB Cross Channel Gain')
