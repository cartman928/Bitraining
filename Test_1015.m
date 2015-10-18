clc
clear


gc(:,1) = [1;1];
gp(:,1) = [1;1];
vc(:,1) = [1.0562 + 0.2200i,0.4991 - 0.4516i];
vp(:,1) = [-0.7415 - 0.2519i;0.3498 + 0.8933i];
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
StepSize = 10^(-5);

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

       

            sum_c1_f(:,1) = H{1,2}*vc(:,2);
            sum_p1_f(:,1) = H{1,2}*vp(:,2);
   
            sum_c2_f(:,1) = [0;0];
            for j = 1:2
                    sum_c2_f(:,1) = sum_c2_f(:,1) + H{1,j}*vc(:,j);
            end

            
            u(:,1) = [0;0];
            for j = 1:2
                    u(:,1) = u(:,1) + H{1,j}*( vc(:,j)*x(iter)+vp(:,j)*xp(iter,j) );
            end
            u(:,1) = u(:,1)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
     

gc(:,1) = gc(:,1)+(StepSize/norm(u(:,1))^2)*u(:,1)*conj(x(iter)-gc(:,1)'*u(:,1))
gp(:,1) = gp(:,1)+(StepSize/norm(u(:,1))^2)*u(:,1)*conj(xp(iter,1)-gp(:,1)'*u(:,1));

end

gc_wiener(:,1) = inv(  H{1,1}*vc(:,1)*vc(:,1)'*H{1,1}' + H{1,1}*vp(:,1)*vp(:,1)'*H{1,1}' + sum_c1_f(:,1)*sum_c1_f(:,1)' +sum_p1_f(:,1)*sum_p1_f(:,1)' + H{1,1}*vc(:,1)*sum_c1_f(:,1)'+sum_c1_f(:,1)*vc(:,1)'*H{1,1}'+eye(2)*sigma^2  ) * ( sum_c2_f(:,1) );    
gp_wiener(:,1) = inv(  H{1,1}*vc(:,1)*vc(:,1)'*H{1,1}' + H{1,1}*vp(:,1)*vp(:,1)'*H{1,1}' + sum_c1_f(:,1)*sum_c1_f(:,1)' +sum_p1_f(:,1)*sum_p1_f(:,1)' + H{1,1}*vc(:,1)*sum_c1_f(:,1)'+sum_c1_f(:,1)*vc(:,1)'*H{1,1}'+eye(2)*sigma^2  ) * ( H{1,1}*vp(:,1) );



