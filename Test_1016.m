%2 user, 2X2 MIMO Channel
%calculate gc(:,1) gp(:,1) gc(:,2) gp(:,2)
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


%{
gc(:,1) = 1;
gp(:,1) = 1;
vc(:,1) = 1;
vp(:,1) = 1;
H{1,1}=-0.9704 + 0.4012i;
%}

sigma = sqrt(10^(-3));
StepSize = 10^(-4);

for iter = 1:2*10^(5) 
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

            
            u(:,k) = [0;0];
            for j = 1:2
                    u(:,k) = u(:,k) + H{k,j}*( vc(:,j)*x(iter)+vp(:,j)*xp(iter,j) );
            end
            u(:,k) = u(:,k)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
     
                  
            gc_wiener(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}' + H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}' + sum_c1_f(:,k)*sum_c1_f(:,k)' +sum_p1_f(:,k)*sum_p1_f(:,k)' + H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2  ) * ( sum_c2_f(:,k) );    
            gp_wiener(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}' + H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}' + sum_c1_f(:,k)*sum_c1_f(:,k)' +sum_p1_f(:,k)*sum_p1_f(:,k)' + H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2  ) * ( H{k,k}*vp(:,k) );

            gc(:,k) = gc(:,k)+StepSize*u(:,k)*conj(x(iter)-gc(:,k)'*u(:,k))
            gp(:,k) = gp(:,k)+StepSize*u(:,k)*conj(xp(iter,k)-gp(:,k)'*u(:,k))

       end

end


%{
gc_wiener_1 =

  -0.3446 - 0.1650i
   0.2692 + 0.6259i


gp_wiener_1 =

   0.1701 + 0.2811i
   0.1678 + 0.1554i

gc_wiener(:,2)


  -0.1469 - 0.0522i
  -0.0879 - 0.3999i

gp_wiener(:,2)


   0.2480 + 0.3436i
  -0.0592 - 0.1199i


gc =

  -0.3350 - 0.1582i  -0.1471 - 0.0523i
   0.2799 + 0.6078i  -0.0884 - 0.4000i


gp =

   0.1750 + 0.2861i   0.2475 + 0.3446i
   0.1774 + 0.1480i  -0.0592 - 0.1174i
%}


