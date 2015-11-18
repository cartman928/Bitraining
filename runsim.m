%% Initialize Parameters

clc
clear

alpha = 0;    %coefficient for block fading model
beta = 0.8^2;  % Attenuation loss from non-direct antennas
n0 = 10^(-2);    %noise variance

Nt = 2;  %Nt antennas for each transmitter
Nr = 2;  %Nr antennas for each receiver
M = 2;   %number of users

upower = ones(1, M);   %power for unicastmpower 
%upower = zeros(1, M);
mpower = ones(1, M);   %power for multicast
%mpower = zeros(1, M);

%upower = sqrt(upower); % Change power to voltage
%mpower = sqrt(mpower); % Change power to voltage

iternums = 1:20; % number of iterations
N_realization = 1000; % Number of times to run simulation

averagerateu = zeros(N_realization, length(iternums));
averageratem = zeros(N_realization, length(iternums));
averagerateu_MaxSINR = zeros(N_realization, length(iternums));
averageratem_MaxSINR = zeros(N_realization, length(iternums));
Eu = zeros(N_realization, length(iternums));
Em = zeros(N_realization, length(iternums));
Eu_MaxSINR = zeros(N_realization, length(iternums));
Em_MaxSINR = zeros(N_realization, length(iternums));

%% Training Length
for traininglength = [10 20] % traininglength 2M
        traininglength
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
    
   
    
    
    %{
    Gu = zeros(Nt, M);
    Gm = zeros(Nt, M);
    Gu_w = zeros(Nt, M);
    Gm_w = zeros(Nt, M);
    Vu = InitialGu; % beamformer unicast
    Vm = InitialGm; % beamformer multicast
    Vu_w = InitialGu; % beamformer unicast
    Vm_w = InitialGm; % beamformer multicast
    %}
    
    
    
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
            %[Vu, Vm] = LS(Z, Gu, Gm, M2, n0, Bbw, BbwBr, upower, mpower); 
            %[Vu_w, Vm_w] = MaxSINR(Z, Gu_w, Gm_w, n0, upower, mpower);
            [Vu, Vm] = LS_cooperation(Z, Gu, Gm, M2, n0, Bbw, BbwBr, upower, mpower);
            [Vu_w, Vm_w] = MaxSINR_cooperation(Z, Gu_w, Gm_w, n0, upower, mpower);
            
            
            %%phase 2: forward training to update receive filter
            [Gu, Gm] = LS(H, Vu, Vm, M1, n0, Bfw, BfwBr, upower, mpower);
            [Gu_w, Gm_w] = MaxSINR(H, Vu_w, Vm_w, n0, upower, mpower);
            
            
            
        averagerateu(realization_idx, numiters, traininglength) = calculate_rateu(H, n0, Vu, Gu, Vm, upower, mpower);
        averageratem(realization_idx, numiters, traininglength) = calculate_ratem(H, n0, Vm, Gm, Vu, upower, mpower);
        averagerateu_MaxSINR(realization_idx, numiters) = calculate_rateu(H, n0, Vu_w, Gu_w, Vm_w, upower, mpower);
        averageratem_MaxSINR(realization_idx, numiters) = calculate_ratem(H, n0, Vm_w, Gm_w, Vu_w, upower, mpower);
        Eu(realization_idx, numiters,traininglength) = MSEu(H, Gu, Gm, Vu, Vm, n0, upower, mpower);
        Em(realization_idx, numiters,traininglength) = MSEm(H, Gu, Gm, Vu, Vm, n0, upower, mpower);
        Eu_MaxSINR(realization_idx, numiters) = MSEu(H, Gu_w, Gm_w, Vu_w, Vm_w, n0, upower, mpower);
        Em_MaxSINR(realization_idx, numiters) = MSEm(H, Gu_w, Gm_w, Vu_w, Vm_w, n0, upower, mpower);
    end
            
    
end

end

%% Plot C(bits/channel)
figure
subplot(2,1,1);
hold on

p1=plot(iternums, mean(averagerateu(:,:,10))+mean(averageratem(:,:,10)),'Color',[0,0.4470,0.7410]);
p2=plot(iternums,mean(averageratem(:,:,10)),'Color',[0,0.4470,0.7410],'Marker','o');
p3=plot(iternums, mean(averagerateu(:,:,10)),'Color',[0,0.4470,0.7410],'Marker','*');

%{
p4=plot(iternums, mean(averagerateu(:,:,16))+mean(averageratem(:,:,16)),'Color',[0.6350,0.0780,0.1840]);
p5=plot(iternums,mean(averageratem(:,:,16)),'Color',[0.6350,0.0780,0.1840],'Marker','o');
p6=plot(iternums, mean(averagerateu(:,:,16)),'Color',[0.6350,0.0780,0.1840],'Marker','*');
%}

p7=plot(iternums, mean(averagerateu(:,:,20))+mean(averageratem(:,:,20)),'Color',[0.8500,0.3250,0.0980]);
p8=plot(iternums,mean(averageratem(:,:,20)),'Color',[0.8500,0.3250,0.0980],'Marker','o');
p9=plot(iternums, mean(averagerateu(:,:,20)),'Color',[0.8500,0.3250,0.0980],'Marker','*');

p10=plot(iternums, mean(averagerateu_MaxSINR)+mean(averageratem_MaxSINR),'k');
p11=plot(iternums, mean(averageratem_MaxSINR),'k','Marker','o');
p12=plot(iternums, mean(averagerateu_MaxSINR),'k','Marker','*');

legend([p10,p1,p7],'C(Max-SINR)',...
                   'C(Bi-Directional Training);2M=10',...    %'C(Bi-Directional Training);2M=16',...
                   'C(Bi-Directional Training);2M=20')
                  
xlabel('Number of iterations')
ylabel('C(bits/channel)')
%title('2 Users;2X2 MIMO Channel;\sigma^2=10^{-2};1000 Realizations;No-Coop')
%title('2 Users;2X2 MIMO Channel;\sigma^2=10^{-2};1000 Realizations;Coop')
title('2 Users;2X2 MIMO Channel;\sigma^2=10^{-2};1000 Realizations')
axis([1 numiters 0 20])

%% Plot MSE
subplot(2,1,2);
hold on

p13=plot(iternums,mean(Em(:,:,10)),'Color',[0,0.4470,0.7410],'Marker','o');
p14=plot(iternums, mean(Eu(:,:,10)),'Color',[0,0.4470,0.7410],'Marker','*');

p15=plot(iternums,mean(Em(:,:,20)),'Color',[0.8500,0.3250,0.0980],'Marker','o');
p16=plot(iternums, mean(Eu(:,:,20)),'Color',[0.8500,0.3250,0.0980],'Marker','*');

p17=plot(iternums, mean(Em_MaxSINR),'k','Marker','o');
p18=plot(iternums, mean(Eu_MaxSINR),'k','Marker','*');

xlabel('Number of iterations')
ylabel('MSE')