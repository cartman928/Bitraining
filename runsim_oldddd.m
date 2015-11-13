%% Initialize Parameters

clc
clear

H(:,:,1,1)=[-0.0482 - 1.2463i   0.4868 - 2.0530i;-0.2119 - 0.1843i  -0.5308 - 0.2194i];
    H(:,:,1,2)=[-0.4955 + 1.3151i  -0.5052 + 0.1263i;-0.4955 + 1.3151i  -0.5052 + 0.1263i];
    H(:,:,2,1)=[0.1879 - 0.2238i   0.6055 + 0.3696i;0.5788 - 0.0921i  -0.2295 - 1.4727i];
    H(:,:,2,2)=[1.0676 + 0.2811i   0.1977 + 0.8887i;1.0676 + 0.2811i   0.1977 + 0.8887i];

N_realization = 100; % Number of times to run simulation 
iter = 200;

%% Start Loop
for realization_idx = 1 : N_realization
        realization_idx

    %% one iteration per block

    InitialGu = [1 1;1 1];    % receive filter unicast
    for k = 1:2
        InitialGu(:, k) = InitialGu(:, k)/norm(InitialGu(:, k));
    end
    
   Gu = InitialGu;
   
    for numiters = 1:iter
        %% bi-directional training
        Bfw = sign(randn(10,2));
        Bbw = sign(randn(10,2));
          
            %%phase 1: backward training to update beamformer
            for sym_idx = 1:10
                for user_idx = 1 : 2
                
                Y1(:,sym_idx,user_idx) = H(:,:,user_idx,user_idx)'*Gu(:,user_idx)*Bbw(sym_idx,user_idx) + sqrt(10^(-3))*(randn(2,1)+1i*randn(2,1))/sqrt(2);
                    for interf_idx = 1 : 2
                            if interf_idx ~= user_idx
                            Y1(:,sym_idx,user_idx) = Y1(:,sym_idx,user_idx) +  H(:,:,interf_idx, user_idx)'*Gu(:,interf_idx)*Bbw(sym_idx,interf_idx);
                            end
                    end
                end 
            end
            
            for user_idx = 1 : 2
            Vu(:,user_idx) = inv(Y1(:,:,user_idx)*Y1(:,:,user_idx)')*Y1(:,:,user_idx)*Bbw(:,user_idx);
            Vu(:,user_idx) = Vu(:,user_idx)/norm(Vu(:,user_idx));
            end
            
            %%phase 2: forward training to update receive filter
            for sym_idx = 1:10
            for user_idx = 1 : 2
                
                Y2(:,sym_idx,user_idx) = H(:,:,user_idx,user_idx)*(Vu(:,user_idx)*Bfw(sym_idx,user_idx)  ) + sqrt(10^(-3))*(randn(2,1)+1i*randn(2,1))/sqrt(2);
                    for interf_idx = 1 : 2
                            if interf_idx ~= user_idx
                            Y2(:,sym_idx,user_idx) = Y2(:,sym_idx,user_idx) +  H(:,:,user_idx,interf_idx)*Vu(:,interf_idx)*Bfw(sym_idx,interf_idx) ;
                            end
                    end
            end
            end
            
            for user_idx = 1 : 2
            Gu(:,user_idx) = inv(Y2(:,:,user_idx)*Y2(:,:,user_idx)')*Y2(:,:,user_idx)*Bfw(:,user_idx);
            Gu(:,user_idx) = Gu(:,user_idx)/norm(Gu(:,user_idx));
            end
            
        %end
        
        
                 for k = 1:2
                 signal = norm(Gu(:,k)'* H(:,:,k,k)*Vu(:,k))^2;
  
                 interference = 0;
                 for interf_idx = 1 : 2
                    if interf_idx ~= k
                        interference = interference + norm(  Gu(:,k)'* H(:,:,k,interf_idx)*(Vu(:,interf_idx))   )^2;
                    end
                 end
    
                 noise = Gu(:,k)'*Gu(:,k)*10^(-3);
                 rate(k) = log2(1+signal/(noise+interference));

                 end
                 averagerate(realization_idx, numiters) = rate(1)+rate(2);

    end
            
    
end

hold on
plot(1:iter, mean(averagerate), 'r');
axis([1 numiters 18.8 21.3])

