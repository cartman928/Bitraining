%2 user, 2X2 MIMO Channel
%Turn Off Common Channel 
%calculate MSE
%LS Filter
%add realization function
clc
clear


sigma = sqrt(10^(-3));

Filterlength = 10; %FilterLength
Realization=100;


for R=1:Realization
    
    R
  
    H(:,:,1,1)=[-0.0482 - 1.2463i   0.4868 - 2.0530i;-0.2119 - 0.1843i  -0.5308 - 0.2194i];
    H(:,:,1,2)=[-0.4955 + 1.3151i  -0.5052 + 0.1263i;-0.4955 + 1.3151i  -0.5052 + 0.1263i];
    H(:,:,2,1)=[0.1879 - 0.2238i   0.6055 + 0.3696i;0.5788 - 0.0921i  -0.2295 - 1.4727i];
    H(:,:,2,2)=[1.0676 + 0.2811i   0.1977 + 0.8887i;1.0676 + 0.2811i   0.1977 + 0.8887i];
        
    for k = 1:2
    gp(:,k)=[1;1];
    gp(:,k)=gp(:,k)/norm(gp(:,k));
    end
    
    InitialGu = [1 1;1 1];    % receive filter unicast
    for k = 1:2
        InitialGu(:, k) = InitialGu(:, k)/norm(InitialGu(:, k));
    end
    
   Gu = InitialGu;


for iteration = 1:200
            %Backward Training
            for iter1 = 1:Filterlength

                    for k =1:2
                        xp_b(k,iter1) = sign(randn);
                    end
                    
                       
                    for k = 1:2
                                neq_p(:,iter1,k) = [0;0];
                                for j = 1:2
                                    if j~=k; 
                                       neq_p(:,iter1,k) = neq_p(:,iter1,k) + H(:,:,j,k)'*gp(:,j)*xp_b(j,iter1);
                                    end
                                end           
                    end
                    
                    nb(:,iter1)=sigma*(1/sqrt(2))*(randn(2,1)+1i*randn(2,1));
                    
                    for k = 1:2    
                                yb(:,iter1,k) = H(:,:,k,k)'*(gp(:,k)*xp_b(k,iter1))...
                                                +neq_p(:,iter1,k)...
                                                +nb(:,iter1);
                    end
                    
            end

                    for k = 1:2
                    vp(:,k)  = inv(yb(:,:,k)*yb(:,:,k)')*yb(:,:,k)*xp_b(k,:)';
                    vp(:,k)=vp(:,k)/norm(vp(:,k));
                    end
            
                    
                    
                    
                    
             %{
              %%phase 1: backward training to update beamformer      
             
                for sym_idx = 1:10
                for user_idx = 1 : 2
                    
                
                Y1(:,sym_idx,user_idx) = H(:,:,user_idx,user_idx)'*Gu(:,user_idx)*xp_b(user_idx,sym_idx) + nb(:,sym_idx);
         
                    for interf_idx = 1 : 2
                            if interf_idx ~= user_idx
                            Y1(:,sym_idx,user_idx) = Y1(:,sym_idx,user_idx) +  H(:,:,interf_idx, user_idx)'*Gu(:,interf_idx)*xp_b(interf_idx,sym_idx);
                            Y1(:,1,1)
                            end
                    end
                end 
                end
                
              
            
           
            for user_idx = 1 : 2
            Vu(:,user_idx) = inv(Y1(:,:,user_idx)*Y1(:,:,user_idx)')*Y1(:,:,user_idx)*xp_b(user_idx,:)';
            Vu(:,user_idx) = Vu(:,user_idx)/norm(Vu(:,user_idx));
            end
                %}
            
           
                    

            %Forward Training  
            for iter2 = 1:Filterlength

                    for k = 1:2
                    xp_f(k,iter2) = sign(randn);  
                    end

                    for k = 1:2
                                neq_p(:,iter2,k) = [0;0];
                                for j = 1:2
                                    if j~=k; 
                                        neq_p(:,iter2,k) = neq_p(:,iter2,k) + H(:,:,k,j)*vp(:,j)*xp_f(j,iter2);
                                    end
                                end
                    end
                    
                    nf(:,iter2)=sigma*(1/sqrt(2))*(randn(2,1)+1i*randn(2,1));
                    
                    for k = 1:2    
                    yf(:,iter2,k) = H(:,:,k,k)*(vp(:,k)*xp_f(k,iter2))...
                                    +neq_p(:,iter2,k)...
                                    +nf(:,iter2);
                    end
                   
            end

                for k = 1:2
                gp(:,k)  = inv(yf(:,:,k)*yf(:,:,k)')*yf(:,:,k)*xp_f(k,:)';
                gp(:,k)=gp(:,k)/norm(gp(:,k));
                end
                
                
                %{
                for sym_idx = 1:10
                for user_idx = 1 : 2
                
                Y2(:,sym_idx,user_idx) = H(:,:,user_idx,user_idx)*(Vu(:,user_idx)*xp_f(user_idx,sym_idx)  ) + nf(:,sym_idx);
                    for interf_idx = 1 : 2
                            if interf_idx ~= user_idx
                            Y2(:,sym_idx,user_idx) = Y2(:,sym_idx,user_idx) +  H(:,:,user_idx,interf_idx)*Vu(:,interf_idx)*xp_f(interf_idx,sym_idx) ;
                            end
                    end
            end
            end
            
            for user_idx = 1 : 2
            Gu(:,user_idx) = inv(Y2(:,:,user_idx)*Y2(:,:,user_idx)')*Y2(:,:,user_idx)*xp_f(user_idx,:)';
            Gu(:,user_idx) = Gu(:,user_idx)/norm(Gu(:,user_idx));
            end
            %}
                

            for k = 1:2
                neq_SINR_p(iteration,k) = 0;
                for j = 1:2
                    if j~=k; 
                    neq_SINR_p(iteration,k) = neq_SINR_p(iteration,k) + norm(gp(:,k)'*H(:,:,k,j)*vp(:,j))^2;
                    end
                end   
            end
            

            for k = 1:2
            SINR_without_stat(iteration,k,R)= ( norm(  gp(:,k)'*H(:,:,k,k)*vp(:,k) )^2   )/...
                                            ( neq_SINR_p(iteration,k) + norm(gp(:,k)'*eye(2)*sigma^2*gp(:,k)) );                                
            end
            C_without_stat(R,iteration)=log2(1+SINR_without_stat(iteration,1,R))+log2(1+SINR_without_stat(iteration,2,R));    
end
           

end
   
hold on
n=1:iteration;

plot(   n,  mean(C_without_stat),'b' )
%legend('C(Bi-Directional Training)')
xlabel('Iteration')
ylabel('C')
title('LS;2 User;2X2 MIMO;Private Messages;Pilot Length 2M=20')
axis([1 iteration 18.8 21.3])


