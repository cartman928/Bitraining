clc
clear


gc(:,1) = [-0.1570 - 1.1223i; 1.9178 + 0.2697i];
gp(:,1) = [ 0.6473 - 1.1667i;-0.3318 + 1.2201i];
vc(:,1) = [1.0562 + 0.2200i,0.4991 - 0.4516i];
vp(:,1) = [-0.7415 - 0.2519i;0.3498 + 0.8933i];
gc(:,2) = [-2.0405 + 0.1186i;0.2665 - 0.1111i];
gp(:,2) = [-0.0690 + 0.2588i;1.4652 + 0.5898i];
vc(:,2) = [-1.3609 + 0.8642i,0.2347 - 0.0695i];
%vc(:,2) = [0;0];
vp(:,2) = [-0.2790 - 0.4776i;0.5612 + 0.9471i];
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

sigma = 10^(-1);
StepSize = 0.5*10^(-5);

for iter = 1:10^(6) 
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

            
            u(:,k) = [0;0];
            for j = 1:2
                    u(:,k) = u(:,k) + H{j,k}.'*( gc(:,j)*x(iter)+gp(:,j)*xp(iter,j) );
            end
            u(:,k) = u(:,k)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];

       end
       
       for k = 1:2
           
           if k == 1
               m=2;
           else
               m=1;
           end 
           
            vc_wiener(:,k) = inv(  H{k,k}.'*gc(:,k)*gc(:,k)'*(H{k,k}').' + H{k,k}.'*gp(:,k)*gp(:,k)'*(H{k,k}').' + sum_c1_b(:,k)*sum_c1_b(:,k)' +sum_p1_b(:,k)*sum_p1_b(:,k)' + H{k,k}.'*gc(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc(:,k)'*(H{k,k}').'+eye(2)*sigma^2  ) * ( sum_c2_b(:,k) );    
            vp_wiener(:,k) = inv(  H{k,k}.'*gc(:,k)*gc(:,k)'*(H{k,k}').' + H{k,k}.'*gp(:,k)*gp(:,k)'*(H{k,k}').' + sum_c1_b(:,k)*sum_c1_b(:,k)' +sum_p1_b(:,k)*sum_p1_b(:,k)' + H{k,k}.'*gc(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc(:,k)'*(H{k,k}').'+eye(2)*sigma^2  ) * ( H{k,k}.'*gp(:,k) );
            vp(:,k) = vp(:,k)+StepSize*u(:,k)*conj(xp(iter,k)-vp(:,k)'*u(:,k));
            
       end
       
       dummy(:,1) = vc(:,1);
       vc(:,1) = vc(:,1)+StepSize*u(:,1)*conj(x(iter)-vc(:,1)'*u(:,1)-vc(:,2)'*u(:,2));
       vc(:,2) = vc(:,2)+StepSize*u(:,2)*conj(x(iter)-vc(:,2)'*u(:,2)-dummy(:,1)'*u(:,1));
       vc

end

%{



%}


