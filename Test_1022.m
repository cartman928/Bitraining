clc
clear


gc(:,1,1) = [1;1];
gp(:,1,1) = [1;1];
vc(:,1,1) = [1.0562 + 0.2200i,0.4991 - 0.4516i];
vp(:,1,1) = [-0.7415 - 0.2519i;0.3498 + 0.8933i];
gc(:,2,1) = [1;1];
gp(:,2,1) = [1;1];
vc(:,2,1) = [-1.3609 + 0.8642i,0.2347 - 0.0695i];
%vc(:,2) = [0;0];
vp(:,2,1) = [-0.2790 - 0.4776i;0.5612 + 0.9471i];
%vp(:,2) = [0;0];
H{1,1}=(1/sqrt(2))*[-0.9704 + 0.4012i (1.9144 - 0.3561i);(0.4516 - 1.4800i) 0.0501 - 0.1627i];
H{1,2}=0.8*(1/sqrt(2))*[-0.6337 + 0.8001i (-1.1485 + 0.2689i);(0.4516 - 1.4800i) -0.6651 - 0.9268i];
H{2,1}=0.8*(1/sqrt(2))*[-0.0165 + 0.9251i (-1.0463 - 0.5763i);(-0.1497 - 1.5829i) -0.7804 - 0.4109i];
H{2,2}=(1/sqrt(2))*[-0.9313 + 0.8060i (0.7313 + 0.0698i);(-0.2850 + 1.1345i) -0.5113 - 0.1662i];


%{
gc(:,1) = 1;
gp(:,1) = 1;
vc(:,1) = 1;
vp(:,1) = 1;
H{1,1}=-0.9704 + 0.4012i;
%}

sigma = 10^(-2);
StepSize = 10^(-4);

for iter = 1:10^(2) 
        if rand-0.5 >= 0
                    x(iter) = 1;
                else
                    x(iter) = -1;
        end

        if rand-0.5 >= 0
                    xp(iter,1) = 1;
                else
                    xp(iter,1) = -1;
        end
        
        if rand-0.5 >= 0
                    xp(iter,2) = 1;
                else
                    xp(iter,2) = -1;
        end

       for k = 1:2 

            sum_c1_f(:,k) = [0;0];
            for j = 1:2
                if j~=k
                    sum_c1_f(:,k) = sum_c1_f(:,k) + H{k,j}*vc(:,j);
                end
            end

            sum_c2_f(:,k) = [0;0];
            for j = 1:2
                    sum_c2_f(:,k) = sum_c2_f(:,k) + H{k,j}*vc(:,j);
            end

            sum_p1_f(:,k) = [0;0];
            for j = 1:2
                if j~=k
                    sum_p1_f(:,k) = sum_p1_f(:,k) + H{k,j}*vp(:,j);
                end
            end

            
             u(:,k,iter) = [0;0];
            for j = 1:2
                    u(:,k,iter) = u(:,k,iter) + H{k,j}*( vc(:,j)*x(iter)+vp(:,j)*xp(iter,j) );
            end
            u(:,k,iter) = u(:,k,iter)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
     
                  
            d_c_estimate(k)=0;
            d_p_estimate(k)=0;
            for n = 1:iter
            d_c_estimate(k) = d_c_estimate(k) + gc(:,k,iter)'*u(:,k,iter-n+1);
            d_p_estimate(k) = d_p_estimate(k) + gp(:,k,iter)'*u(:,k,iter-n+1);
            end
            
            for n = 1:iter
            gc(:,k,n+1) = gc(:,k,n)+StepSize*u(:,k,iter-n+1)*conj(x(iter)-d_c_estimate(k))
            gp(:,k,n+1) = gp(:,k,n)+StepSize*u(:,k,iter-n+1)*conj(xp(iter,k)-d_p_estimate(k));
            end

       end

end



for iter = 1:10^(2) 
        if rand-0.5 >= 0
                    x(iter) = 1;
                else
                    x(iter) = -1;
        end

        if rand-0.5 >= 0
                    xp(iter,1) = 1;
                else
                    xp(iter,1) = -1;
        end
        
        if rand-0.5 >= 0
                    xp(iter,2) = 1;
                else
                    xp(iter,2) = -1;
        end

       for k = 1:2 

            sum_c1_b(:,k) = [0;0];
            for j = 1:2
                if j~=k
                    sum_c1_b(:,k) = sum_c1_b(:,k) + H{j,k}.'*gc(:,j);
                end
            end

            sum_c2_b(:,k) = [0;0];
            for j = 1:2
                    sum_c2_b(:,k) = sum_c2_b(:,k) + H{j,k}.'*gc(:,j);
            end

            sum_p1_b(:,k) = [0;0];
            for j = 1:2
                if j~=k
                    sum_p1_b(:,k) = sum_p1_b(:,k) + H{j,k}.'*gp(:,j);
                end
            end

            
            u(:,k,iter) = [0;0];
            for j = 1:2
                    u(:,k,iter) = u(:,k) + H{j,k}.'*( gc(:,j,100)*x(iter)+gp(:,j,100)*xp(iter,j) );
            end
            u(:,k,iter) = u(:,k,iter)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
     
                  
            vc_wiener(:,k) = inv(  H{k,k}.'*gc(:,k)*gc(:,k)'*(H{k,k}').' + H{k,k}.'*gp(:,k)*gp(:,k)'*(H{k,k}').' + sum_c1_b(:,k)*sum_c1_b(:,k)' +sum_p1_b(:,k)*sum_p1_b(:,k)' + H{k,k}.'*gc(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc(:,k)'*(H{k,k}').'+eye(2)*sigma^2  ) * ( sum_c2_b(:,k) );    
            vp_wiener(:,k) = inv(  H{k,k}.'*gc(:,k)*gc(:,k)'*(H{k,k}').' + H{k,k}.'*gp(:,k)*gp(:,k)'*(H{k,k}').' + sum_c1_b(:,k)*sum_c1_b(:,k)' +sum_p1_b(:,k)*sum_p1_b(:,k)' + H{k,k}.'*gc(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc(:,k)'*(H{k,k}').'+eye(2)*sigma^2  ) * ( H{k,k}.'*gp(:,k) );

            
            d_c_estimate(k)=0;
            d_p_estimate(k)=0;
            for n = 1:iter
            d_c_estimate(k) = d_c_estimate(k) + vc(:,k,iter)'*u(:,k,iter-n+1);
            d_p_estimate(k) = d_p_estimate(k) + vp(:,k,iter)'*u(:,k,iter-n+1);
            end
            
            for n = 1:iter
            vc(:,k,n+1) = vc(:,k,n)+StepSize*u(:,k,iter-n+1)*conj(x(iter)-d_c_estimate(k))
            vp(:,k,n+1) = vp(:,k,n)+StepSize*u(:,k,iter-n+1)*conj(xp(iter,k)-d_p_estimate(k));
            end
            
            

       end

end

SINR_C(1)=(   norm( gc(:,1)'*H{1,1}*vc(:,1)+gc(:,2)'*H{1,2}*vc(:,2) )^2   )/( norm( gc(:,1)'*sigma^2*gc(:,1) ) +  norm( gc(:,1)'*H{1,1}*vp(:,1)+gc(:,1)'*H{1,2}*vp(:,2) )^2  );
SINR_C(2)=(   norm( gc(:,2)'*H{2,1}*vc(:,1)+gc(:,2)'*H{2,2}*vc(:,2) )^2   )/( norm( gc(:,2)'*sigma^2*gc(:,2) ) +  norm( gc(:,2)'*H{2,1}*vp(:,1)+gc(:,2)'*H{2,2}*vp(:,2) )^2  );

SINR_P(2)=(   norm( gp(:,1)'*H{1,1}*vp(:,1) )^2   )/( norm( gp(:,1)'*sigma^2*gp(:,1) ) +  norm( gp(:,1)'*H{1,1}*vc(:,1)+gp(:,1)'*H{1,2}*vc(:,2)+gp(:,1)'*H{1,2}*vp(:,2) )^2  );
SINR_P(2)=(   norm( gp(:,2)'*H{2,2}*vp(:,2) )^2   )/( norm( gp(:,2)'*sigma^2*gp(:,2) ) +  norm( gp(:,2)'*H{2,1}*vc(:,1)+gp(:,2)'*H{2,2}*vc(:,2)+gp(:,2)'*H{2,1}*vp(:,1) )^2  );

C=0;
for i = 1:2
    
    C=log2(1+SINR_C(i))+log2(1+SINR_P(i))
    
end