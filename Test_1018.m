%2 user, 2X2 MIMO Channel
%calculate vc(:,1) vp(:,1) vc(:,2) vp(:,2) with cooperation
clc
clear


gc(:,1,1) = [-0.3520 - 0.1639i;0.2730 + 0.6405i];
gp(:,1,1) = [0.1697 + 0.2841i;0.1733 + 0.1567i];
gc(:,2,1) = [-0.1478 - 0.0519i;-0.0879 - 0.4007i];
gp(:,2,1) = [0.2490 + 0.3453i;-0.0598 - 0.1202i];

gc(:,1)=gc(:,1)/norm(gc(:,1));
gc(:,2)=gc(:,2)/norm(gc(:,2));
gp(:,1)=gp(:,1)/norm(gp(:,1));
gp(:,2)=gp(:,2)/norm(gp(:,2));


vc(:,1) = [1.0562 + 0.2200i,0.4991 - 0.4516i];
vp(:,1) = [-0.7415 - 0.2519i;0.3498 + 0.8933i];
vc(:,2) = [-1.3609 + 0.8642i,0.2347 - 0.0695i];
vp(:,2) = [-0.2790 - 0.4776i;0.5612 + 0.9471i];

H{1,1}=(1/sqrt(2))*[-0.9704 + 0.4012i (1.9144 - 0.3561i);(0.4516 - 1.4800i) 0.0501 - 0.1627i];
H{1,2}=0.8*(1/sqrt(2))*[-0.6337 + 0.8001i (-1.1485 + 0.2689i);(0.4516 - 1.4800i) -0.6651 - 0.9268i];
H{2,1}=0.8*(1/sqrt(2))*[-0.0165 + 0.9251i (-1.0463 - 0.5763i);(-0.1497 - 1.5829i) -0.7804 - 0.4109i];
H{2,2}=(1/sqrt(2))*[-0.9313 + 0.8060i (0.7313 + 0.0698i);(-0.2850 + 1.1345i) -0.5113 - 0.1662i];



sigma = sqrt(10^(-3));
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

vc_wiener =

   0.0210 - 0.2243i   0.1323 - 0.0615i
   0.2431 - 0.0885i  -0.3279 - 0.0058i

vp_wiener =

  -0.0173 + 0.1862i  -0.1159 + 0.1133i
   0.0264 - 0.3804i  -0.3255 - 0.1640i

vc =

   0.0938 - 0.1661i   0.1053 - 0.0461i
   0.1833 + 0.0608i  -0.0980 + 0.2426i

vp =

  -0.0171 + 0.1863i  -0.1158 + 0.1134i
   0.0262 - 0.3801i  -0.3250 - 0.1631i
%}


