%2 user, 2X2 MIMO Channel
%Turn Off Common Channel 
%calculate MSE
%LS Filter
%add realization function
clc
clear



sigma = sqrt(10^(-3));

for i = [4,10]; %FilterLength
    
    i

for k=1:3
SINR_without_stat(:,k,i)= zeros(100,1);
SINR_know_stat(:,k,i)= zeros(100,1);
end

Realization=100;


 for R=1:Realization
     
        R
        
        H{1,1}=(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
        H{1,2}=0.8*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
        H{2,1}=0.8*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
        H{2,2}=(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
        H{3,3}=(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
        H{1,3}=0.8*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
        H{2,3}=0.8*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
        H{3,1}=0.8*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
        H{3,2}=0.8*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];

        Z{1,1}=H{1,1}.';
        Z{1,2}=H{2,1}.';
        Z{2,1}=H{1,2}.';
        Z{2,2}=H{2,2}.';
        Z{3,3}=H{3,3}.';
        Z{3,1}=H{1,3}.';
        Z{3,2}=H{2,3}.';
        Z{1,3}=H{3,1}.';
        Z{2,3}=H{3,2}.';
        
        for k = 1:3
        gp(:,k)=[1;1];
        gp_w(:,k)=[1;1];
        gp(:,k)=gp(:,k)/norm(gp(:,k));
        gp_w(:,k)=gp_w(:,k)/norm(gp_w(:,k));
        end


for iteration = 1:50

    
            for k = 1:3
            gp(:,k)=gp(:,k)/norm(gp(:,k));
            gp_w(:,k)=gp_w(:,k)/norm(gp_w(:,k));
            end
            
           
            %Backward Training
            for iter1 = 1:i

                    for k =1:3
                        if rand-0.5 >= 0
                                    xp_b(k,iter1) = 1;
                                else
                                    xp_b(k,iter1) = -1;
                        end
                    end
                    
                       
                    for k = 1:3

                                neq_p(:,iter1,k) = [0;0];
                                for j = 1:2
                                    if j~=k; 
                                        neq_p(:,iter1,k) = neq_p(:,iter1,k) + Z{k,j}*gp(:,j)*xp_b(j,iter1);
                                    end
                                end

                                neq_p_w(:,:,k) = [0 0;0 0];
                                for j = 1:2
                                    if j~=k;
                                        neq_p_w(:,:,k) = neq_p_w(:,:,k) + Z{k,j}*gp_w(:,j)*gp_w(:,j)'*Z{k,j}';
                                    end
                                end
                                
                    end

                                  
                    for k = 1:3    
                                yb(:,iter1,k) = Z{k,k}*(gp(:,k)*xp_b(k,iter1))...
                                                +neq_p(:,iter1,k)...
                                                +sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))];
                    end
                    
            end
           
                    for k = 1:3
                            vp_w(:,k) = inv(  Z{k,k}*gp_w(:,k)*gp_w(:,k)'*Z{k,k}'... 
                                            + neq_p_w(:,:,k)*(neq_p_w(:,:,k))'...  
                                            + eye(2)*sigma^2  ) * (Z{k,k}*gp_w(:,k)); 
                    end
            
                    for k = 1:3
                    vp(:,k)  = inv(yb(:,:,k)*yb(:,:,k)')*yb(:,:,k)*xp_b(k,:)';
                    end
            
                    for k = 1:3
                    vp(:,k)=vp(:,k)/norm(vp(:,k));
                    vp_w(:,k)=vp_w(:,k)/norm(vp_w(:,k));
                    end
                    
           
       

            %Forward Training  
            for iter2 = 1:i

                    for k = 1:3
                        if rand-0.5 >= 0
                                    xp_f(k,iter2) = 1;
                                else
                                    xp_f(k,iter2) = -1;
                        end
                    end

                    for k = 1:3

                                neq_p(:,iter2,k) = [0;0];
                                for j = 1:2
                                    if j~=k; 
                                        neq_p(:,iter2,k) = neq_p(:,iter2,k) + H{k,j}*vp(:,j)*xp_f(j,iter2);
                                    end
                                end

                                neq_p_w(:,:,k) = [0 0;0 0];
                                for j = 1:2
                                    if j~=k;
                                        neq_p_w(:,:,k) = neq_p_w(:,:,k) + H{k,j}*vp_w(:,j)*vp_w(:,j)'*H{k,j}';
                                    end
                                end

                    end
                    
                    for k = 1:3    
                    yf(:,iter2,k) = H{k,k}*(vp(:,k)*xp_f(k,iter2))...
                                    +neq_p(:,iter2,k)...
                                    +sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))];
                    end
                   
            end
            
    
                for k = 1:3 
                gp_w(:,k) = inv(  H{k,k}*vp_w(:,k)*vp_w(:,k)'*H{k,k}'...
                                + neq_p_w(:,:,k)*(neq_p_w(:,:,k))'...  
                                + eye(2)*sigma^2  ) * (H{k,k}*vp_w(:,k)); 
                end

                for k = 1:3
                gp(:,k)  = inv(yf(:,:,k)*yf(:,:,k)')*yf(:,:,k)*xp_f(k,:)';
                end
                
                      
           
    
            for k = 1:3
        
                %%%%%%%%%%%
                neq_SINR_p(iteration,k) = 0;
                for j = 1:3
                    if j~=k; 
                neq_SINR_p(iteration,k) = neq_SINR_p(iteration,k) + norm(gp(:,k)'*H{k,j}*vp(:,j))^2;
                    end
                end
                %%%%%%%%%%%
                neq_SINR_p_w(iteration,k) = 0;
                for j = 1:3
                    if j~=k; 
                neq_SINR_p_w(iteration,k) = neq_SINR_p_w(iteration,k) + norm(gp_w(:,k)'*H{k,j}*vp_w(:,j))^2;
                    end
                end
        
            end
    
            for k = 1:3
            SINR_without_stat(iteration,k,i)= SINR_without_stat(iteration,k,i)+(norm(  gp(:,k)'*H{k,k}*vp(:,k) )^2/( neq_SINR_p(iteration,k) + norm( gp(:,k)'*eye(2)*sigma^2*gp(:,k) )    ))/Realization;
            SINR_know_stat(iteration,k,i)= SINR_know_stat(iteration,k,i)+(norm(   gp_w(:,k)'*H{k,k}*vp_w(:,k) )^2/( neq_SINR_p_w(iteration,k) + norm( gp_w(:,k)'*eye(2)*sigma^2*gp_w(:,k) )   ))/Realization;
            end
           
    
    end
           
 end
 
 n=1:iteration;
 plot(   n,  log2(1+SINR_without_stat(n,1,i))+log2(1+SINR_without_stat(n,2,i))+log2(1+SINR_without_stat(n,3,i)) );
 hold on
 
 end

plot( n,log2(1+SINR_know_stat(n,1,i))+log2(1+SINR_know_stat(n,2,i))+log2(1+SINR_know_stat(n,3,i)) );

legend('C(Bi-Directional);2M=4',...
       'C(Bi-Directional);2M=8','C(Bi-Directional);2M=16',...
       'C(Max-SINR)')
xlabel('Iteration')
ylabel('C')
title('LS;3 User;2X2 MIMO;Private Messages;1000 Realization')
axis([1 iteration 0 25])


