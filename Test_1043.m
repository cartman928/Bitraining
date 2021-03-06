%3 user, 2X2 MIMO Channel
%All Channels(without cooperation)
%calculate MSE
%LS Filter
%add realization function
clc
clear


sigma = sqrt(10^(-3));

FilterLength = 10; %FilterLength
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

        Z{1,1}=H{1,1}';
        Z{1,2}=H{2,1}';
        Z{2,1}=H{1,2}';
        Z{2,2}=H{2,2}';
        Z{3,3}=H{3,3}';
        Z{3,1}=H{1,3}';
        Z{3,2}=H{2,3}';
        Z{1,3}=H{3,1}';
        Z{2,3}=H{3,2}';
        
  
        
        

        for k = 1:3
            gp(:,k)=[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
            gp_w(:,k)=[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
            gc(:,k)=[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
            gc_w(:,k)=[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
            gc(:,k)=gc(:,k)/norm(gc(:,k));
            gc_w(:,k)=gc_w(:,k)/norm(gc_w(:,k));
            gp(:,k)=gp(:,k)/norm(gp(:,k));
            gp_w(:,k)=gp_w(:,k)/norm(gp_w(:,k));
        end


        for iteration = 1:20
    
           
            %Backward Training
            for iter1 = 1:FilterLength

                    for k =1:3
                        if rand-0.5 >= 0
                                    xp_b(k,iter1) = 1;
                                else
                                    xp_b(k,iter1) = -1;
                        end
                    end
                    
                    if rand-0.5 >= 0
                                xc_b(iter1) = 1;
                            else
                                xc_b(iter1) = -1;
                    end
                    
                       
                    for k = 1:3

                                neq_p(:,iter1,k) = [0;0];
                                for j = 1:3
                                    if j~=k; 
                                        neq_p(:,iter1,k) = neq_p(:,iter1,k) + Z{k,j}*gp(:,j)*xp_b(j,iter1);
                                    end
                                end

                                neq_p_w(:,:,k) = [0 0;0 0];
                                for j = 1:3
                                    if j~=k;
                                        neq_p_w(:,:,k) = neq_p_w(:,:,k) + Z{k,j}*gp_w(:,j)*gp_w(:,j)'*Z{k,j}';
                                    end
                                end
                                
                                neq(:,iter1,k) = [0;0];
                                for j = 1:3
                                    if j~=k; 
                                        neq(:,iter1,k) = neq(:,iter1,k) + Z{k,j}*gc(:,j)*xc_b(iter1);
                                    end
                                end

                                neq_w(:,iter1,k) = [0;0];
                                for j = 1:3
                                    if j~=k;
                                        neq_w(:,iter1,k) = neq_w(:,iter1,k) + Z{k,j}*gc_w(:,j);
                                    end
                                end

                                all_w(:,iter1,k) = [0;0];
                                for j = 1:3
                                        all_w(:,iter1,k) = all_w(:,iter1,k) + Z{k,j}*gc_w(:,j);
                                end

                           
                    end
                    
                    for k = 1:3 
                                yb(:,iter1,k) = Z{k,k}*(gp(:,k)*xp_b(k,iter1))...
                                                +neq_p(:,iter1,k)...
                                                +sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))]...
                                                +Z{k,k}*(gc(:,k)*xc_b(iter1))...
                                                +neq(:,iter1,k);
                    end
                           
                    
                    
            end
            
                    for k = 1:3
                            vp_w(:,k) = inv(  Z{k,k}*gp_w(:,k)*gp_w(:,k)'*Z{k,k}'... 
                                            + neq_p_w(:,:,k)...  
                                            + eye(2)*sigma^2 ...
                                            + Z{k,k}*gc_w(:,k)*gc_w(:,k)'*Z{k,k}'...
                                            + neq_w(:,iter1,k)*(neq_w(:,iter1,k))'...
                                            + Z{k,k}*gc_w(:,k)*(neq_w(:,iter1,k))'...
                                            + (Z{k,k}*gc_w(:,k)*(neq_w(:,iter1,k))')'...
                                            ) * (Z{k,k}*gp_w(:,k)); 
                                        
                            vc_w(:,k) = inv(  Z{k,k}*gc_w(:,k)*gc_w(:,k)'*Z{k,k}'  + neq_w(:,iter1,k)*(neq_w(:,iter1,k))'...  
                                        + Z{k,k}*gc_w(:,k)*(neq_w(:,iter1,k))'+ (Z{k,k}*gc_w(:,k)*(neq_w(:,iter1,k))')'...
                                        + eye(2)*sigma^2 ...
                                        + Z{k,k}*gp_w(:,k)*gp_w(:,k)'*Z{k,k}'...
                                        + neq_p_w(:,:,k)...
                                        ) * all_w(:,iter1,k); 
                    end
            
                    for k = 1:3
                    vp(:,k)  = inv(yb(:,:,k)*yb(:,:,k)')*yb(:,:,k)*xp_b(k,:)';
                    vc(:,k)  = inv(yb(:,:,k)*yb(:,:,k)')*yb(:,:,k)*xc_b';
                    end
            
                    for k = 1:3
                    vp(:,k)=vp(:,k)/norm(vp(:,k));
                    vp_w(:,k)=vp_w(:,k)/norm(vp_w(:,k));
                    vc(:,k)=vc(:,k)/norm(vc(:,k));
                    vc_w(:,k)=vc_w(:,k)/norm(vc_w(:,k));
                    end
                    
       

            %Forward Training  
            for iter2 = 1:FilterLength

                    for k = 1:3
                        if rand-0.5 >= 0
                                    xp_f(k,iter2) = 1;
                                else
                                    xp_f(k,iter2) = -1;
                        end
                    end
                    
                    if rand-0.5 >= 0
                                xc_f(iter2) = 1;
                            else
                                xc_f(iter2) = -1;
                    end

                    for k = 1:3

                                neq_p(:,iter2,k) = [0;0];
                                for j = 1:3
                                    if j~=k; 
                                        neq_p(:,iter2,k) = neq_p(:,iter2,k) + H{k,j}*vp(:,j)*xp_f(j,iter2);
                                    end
                                end

                                neq_p_w(:,:,k) = [0 0;0 0];
                                for j = 1:3
                                    if j~=k;
                                        neq_p_w(:,:,k) = neq_p_w(:,:,k) + H{k,j}*vp_w(:,j)*vp_w(:,j)'*H{k,j}';
                                    end
                                end
                                
                                neq(:,iter2,k) = [0;0];
                                for j = 1:3
                                    if j~=k; 
                                        neq(:,iter2,k) = neq(:,iter2,k) + H{k,j}*vc(:,j)*xc_f(iter2);
                                    end
                                end

                                neq_w(:,iter2,k) = [0;0];
                                for j = 1:3
                                    if j~=k;
                                        neq_w(:,iter2,k) = neq_w(:,iter2,k) + H{k,j}*vc_w(:,j);
                                    end
                                end

                                all_w(:,iter2,k) = [0;0];
                                for j = 1:3
                                        all_w(:,iter2,k) = all_w(:,iter2,k) + H{k,j}*vc_w(:,j);
                                end

                    end
                    
                    for k = 1:3  
                    yf(:,iter2,k) = H{k,k}*(vp(:,k)*xp_f(k,iter2))...
                                    +neq_p(:,iter2,k)...
                                    +sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))]...
                                    +H{k,k}*(vc(:,k)*xc_f(iter2))...
                                    +neq(:,iter2,k);
                    end
                   
            end 
    
                for k = 1:3
                gp_w(:,k) = inv(  H{k,k}*vp_w(:,k)*vp_w(:,k)'*H{k,k}'...
                                + neq_p_w(:,:,k)...  
                                + eye(2)*sigma^2 ...
                                + H{k,k}*vc_w(:,k)*vc_w(:,k)'*H{k,k}'...
                                + neq_w(:,iter2,k)*(neq_w(:,iter2,k))'...
                                + neq_w(:,iter2,k)*vc_w(:,k)'*H{k,k}'...
                                + H{k,k}*vc_w(:,k)*(neq_w(:,iter2,k))'...
                                ) * (H{k,k}*vp_w(:,k));
                            
                gc_w(:,k) = inv(  H{k,k}*vc_w(:,k)*vc_w(:,k)'*H{k,k}'  + neq_w(:,iter2,k)*(neq_w(:,iter2,k))'...  
                            + neq_w(:,iter2,k)*vc_w(:,k)'*H{k,k}'  + H{k,k}*vc_w(:,k)*(neq_w(:,iter2,k))'...
                            + eye(2)*sigma^2 ...
                            + H{k,k}*vp_w(:,k)*vp_w(:,k)'*H{k,k}'...
                            + neq_p_w(:,:,k)...
                            ) * all_w(:,iter2,k);  
                        
          
                
                            
                end

                for k = 1:3
                gp(:,k)  = inv(yf(:,:,k)*yf(:,:,k)')*yf(:,:,k)*xp_f(k,:)';
                gc(:,k)  = inv(yf(:,:,k)*yf(:,:,k)')*yf(:,:,k)*xc_f';
               
                
                end
                
                for k = 1:3
                gp(:,k)=gp(:,k)/norm(gp(:,k));
                gp_w(:,k)=gp_w(:,k)/norm(gp_w(:,k));
                gc(:,k)=gc(:,k)/norm(gc(:,k));
                gc_w(:,k)=gc_w(:,k)/norm(gc_w(:,k));
                end      
                
            
    
            %for SINR
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
                
                all_SINR(:,iteration,k) = [0;0];
                for j = 1:3
                all_SINR(:,iteration,k) = all_SINR(:,iteration,k) + H{k,j}*vc(:,j);
                end

                all_SINR_w(:,iteration,k) = [0;0];
                for j = 1:3
                all_SINR_w(:,iteration,k) = all_SINR_w(:,iteration,k) + H{k,j}*vc_w(:,j);
                end
                
                %%%%%%%%%%
                all_SINR_cp(iteration,k) = 0;
                for j = 1:3
                all_SINR_cp(iteration,k) = all_SINR_cp(iteration,k) + norm(gc(:,k)'*H{k,j}*vp(:,j))^2;
                end
                %%%%%%%%%%
                all_SINR_w_cp(iteration,k) = 0;
                for j = 1:3
                all_SINR_w_cp(iteration,k) = all_SINR_w_cp(iteration,k) + norm(gc_w(:,k)'*H{k,j}*vp_w(:,j))^2;
                end
                
                all_SINR_pc(iteration,k) = 0;
                for j = 1:3
                all_SINR_pc(iteration,k) = all_SINR_pc(iteration,k) + gp(:,k)'*H{k,j}*vc(:,j);
                end
                
                all_SINR_w_pc(iteration,k) = 0;
                for j = 1:3
                all_SINR_w_pc(iteration,k) = all_SINR_w_pc(iteration,k) + gp_w(:,k)'*H{k,j}*vc_w(:,j);
                end
        
            end
    
            for k = 1:3
            SINR_p_without_stat(iteration,k,R)= ...
                      (    norm(  gp(:,k)'*H{k,k}*vp(:,k) )^2/(  neq_SINR_p(iteration,k)+ norm( gp(:,k)'*eye(2)*sigma^2*gp(:,k) )  + norm(all_SINR_pc(iteration,k))^2    ));
            SINR_p_know_stat(iteration,k,R)= ...
                      (    norm(  gp_w(:,k)'*H{k,k}*vp_w(:,k) )^2/(  neq_SINR_p_w(iteration,k) + norm( gp_w(:,k)'*eye(2)*sigma^2*gp_w(:,k) ) + norm(all_SINR_w_pc(iteration,k))^2  ));    
            SINR_c_without_stat(iteration,k,R)= ...
                      (    norm(gc(:,k)'*all_SINR(:,iteration,k)   )^2/(  norm( gc(:,k)'*eye(2)*sigma^2*gc(:,k) )  + all_SINR_cp(iteration,k)     ));
            SINR_c_know_stat(iteration,k,R)= ...
                      (    norm( gc_w(:,k)'*all_SINR_w(:,iteration,k)  )^2/(  norm( gc_w(:,k)'*eye(2)*sigma^2*gc_w(:,k) ) + all_SINR_w_cp(iteration,k)    )); 
            end
            
            
            C_c_without_stat(R,iteration)=0;
            C_c_know_stat(R,iteration)=0;
            C_p_without_stat(R,iteration)=0;
            C_p_know_stat(R,iteration)=0;
            for k = 1:3
            C_c_without_stat(R,iteration)=C_c_without_stat(R,iteration)+log2(1+SINR_c_without_stat(iteration,k,R));
            C_c_know_stat(R,iteration)=C_c_know_stat(R,iteration)+log2(1+SINR_c_know_stat(iteration,k,R));
            C_p_without_stat(R,iteration)=C_p_without_stat(R,iteration)+log2(1+SINR_p_without_stat(iteration,k,R));
            C_p_know_stat(R,iteration)=C_p_know_stat(R,iteration)+log2(1+SINR_p_know_stat(iteration,k,R));
            end
 
    end
           
end
   

n=1:iteration;

plot(   n,  mean(C_c_without_stat)+mean(C_p_without_stat),   'r', n,  mean(C_c_know_stat)+mean(C_p_know_stat), 'y' )


legend('C(Bi-Directional Training)','C(Max-SINR)')
xlabel('Iteration')
ylabel('C')
title('LS;3 User;Fixed 2X2 MIMO;Pilot Length 2M=20;no coop')
axis([1 iteration 10 17])
