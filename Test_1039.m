%2 user, 2X2 MIMO Channel
%Turn Off Private Channel 
%calculate MSE
%LS Filter
clc
clear

H{1,1}=(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
H{1,2}=0.8*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
H{2,1}=0.8*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];
H{2,2}=(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];

Z{1,1}=H{1,1}.';
Z{1,2}=H{2,1}.';
Z{2,1}=H{1,2}.';
Z{2,2}=H{2,2}.';

sigma = sqrt(10^(-3));
StepSize = 10^(-3);

i = 40; %FilterLength

for iteration = 1:20

    iteration
    
    gc(:,1)=[1;1];
    gc(:,2)=[1;1];
    gc_w(:,1)=[1;1];
    gc_w(:,2)=[1;1];
    gc(:,1)=gc(:,1)/norm(gc(:,1));
    gc(:,2)=gc(:,2)/norm(gc(:,2));
    gc_w(:,1)=gc_w(:,1)/norm(gc_w(:,1));
    gc_w(:,2)=gc_w(:,2)/norm(gc_w(:,2));
    
    for loop=1:iteration
    
            
            for k = 1:2
            gc(:,k)=gc(:,k)/norm(gc(:,k));
            gc_w(:,k)=gc_w(:,k)/norm(gc_w(:,k));
            end
            
           
            %Backward Training
            for iter1 = 1:i

                    if rand-0.5 >= 0
                                xc_b(iter1) = 1;
                            else
                                xc_b(iter1) = -1;
                    end
                    
                       
                    for k = 1:2

                                neq(:,iter1,k) = [0;0];
                                for j = 1:2
                                    if j~=k; 
                                        neq(:,iter1,k) = neq(:,iter1,k) + Z{k,j}*gc(:,j)*xc_b(iter1);
                                    end
                                end

                                neq_w(:,iter1,k) = [0;0];
                                for j = 1:2
                                    if j~=k;
                                        neq_w(:,iter1,k) = neq_w(:,iter1,k) + Z{k,j}*gc_w(:,j);
                                    end
                                end

                                all_w(:,iter1,k) = [0;0];
                                for j = 1:2
                                        all_w(:,iter1,k) = all_w(:,iter1,k) + Z{k,j}*gc_w(:,j);
                                end

                    end
                    
                    for k = 1:2    
                                yb(:,iter1,k) = Z{k,k}*(gc(:,k)*xc_b(iter1))...
                                                +neq(:,iter1,k)...
                                                +sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))];
                    end
                            
                    
            end
            
            for k = 1:2 
                    vc_w(:,k) = inv(  Z{k,k}*gc_w(:,k)*gc_w(:,k)'*Z{k,k}'  + neq_w(:,iter1,k)*(neq_w(:,iter1,k))'...  
                                    + Z{k,k}*gc_w(:,k)*(neq_w(:,iter1,k))'+ (Z{k,k}*gc_w(:,k)*(neq_w(:,iter1,k))')'...
                                    + eye(2)*sigma^2  ) * all_w(:,iter1,k); 
            end
            
            for k = 1:2
            vc(:,k)  = inv(yb(:,:,k)*yb(:,:,k)')*yb(:,:,k)*xc_b';
            end
            
            for k = 1:2
            vc(:,k)=vc(:,k)/norm(vc(:,k));
            vc_w(:,k)=vc_w(:,k)/norm(vc_w(:,k));
            end
    
                     

            %Forward Training  
            for iter2 = 1:i

                    if rand-0.5 >= 0
                                xc_f(iter2) = 1;
                            else
                                xc_f(iter2) = -1;
                    end

                    for k = 1:2

                                neq(:,iter2,k) = [0;0];
                                for j = 1:2
                                    if j~=k; 
                                        neq(:,iter2,k) = neq(:,iter2,k) + H{k,j}*vc(:,j)*xc_f(iter2);
                                    end
                                end

                                neq_w(:,iter2,k) = [0;0];
                                for j = 1:2
                                    if j~=k;
                                        neq_w(:,iter2,k) = neq_w(:,iter2,k) + H{k,j}*vc_w(:,j);
                                    end
                                end

                                all_w(:,iter2,k) = [0;0];
                                for j = 1:2
                                        all_w(:,iter2,k) = all_w(:,iter2,k) + H{k,j}*vc_w(:,j);
                                end

                    end
                    
                    for k = 1:2    
                    yf(:,iter2,k) = H{k,k}*(vc(:,k)*xc_f(iter2))...
                                    +neq(:,iter2,k)...
                                    +sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))];
                    end
                   
            end 
    
            for k = 1:2 
            gc_w(:,k) = inv(  H{k,k}*vc_w(:,k)*vc_w(:,k)'*H{k,k}'  + neq_w(:,iter2,k)*(neq_w(:,iter2,k))'...  
                            + neq_w(:,iter2,k)*vc_w(:,k)'*H{k,k}'  + H{k,k}*vc_w(:,k)*(neq_w(:,iter2,k))'...
                            + eye(2)*sigma^2  ) * all_w(:,iter2,k); 
            end
            
            for k = 1:2
            gc(:,k)  = inv(yf(:,:,k)*yf(:,:,k)')*yf(:,:,k)*xc_f';
            end
                      
    end
    
    for k = 1:2
        
        all_SINR(:,iteration,k) = [0;0];
        for j = 1:2
        all_SINR(:,iteration,k) = all_SINR(:,iteration,k) + H{k,j}*vc(:,j);
        end

        all_SINR_w(:,iteration,k) = [0;0];
        for j = 1:2
        all_SINR_w(:,iteration,k) = all_SINR_w(:,iteration,k) + H{k,j}*vc_w(:,j);
        end
        
    end
    
    for k = 1:2
    SINR_without_stat(iteration,k)= norm(gc(:,k)'*all_SINR(:,iteration,k))^2/norm( gc(:,k)'*eye(2)*sigma^2*gc(:,k) );
    SINR_know_stat(iteration,k)= norm(( gc_w(:,k)'*all_SINR_w(:,iteration,k)))^2/norm( gc_w(:,k)'*eye(2)*sigma^2*gc_w(:,k) ); 
    end
    
           
end
   

n=1:iteration;

plot(   n,log2(1+SINR_without_stat(n,1))+log2(1+SINR_without_stat(n,2)),n,log2(1+SINR_know_stat(n,1))+log2(1+SINR_know_stat(n,2)))
legend('C(Bi-Directional Training)','C(Max-SINR)')
xlabel('Iteration')
ylabel('C')
title('LS;1 User;Fixed 2X2 MIMO;Pilot Length=50')
axis([1 iteration 0 30])


