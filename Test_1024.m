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

vc(:,1) = [0.3902 + 0.2111i;0.7619 - 0.2572i];
vp(:,1) = [-0.3188 + 0.6144i;0.6025 + 0.8855i];
vc(:,2) = [0.4874 - 0.1169i;0.5935 + 0.2081i];
vp(:,2) = [-0.1296 - 0.1773i;0.3754 + 0.5062i];

vc(:,1)=vc(:,1)/norm(vc(:,1));
vc(:,2)=vc(:,2)/norm(vc(:,2));
vp(:,1)=vp(:,1)/norm(vp(:,1));
vp(:,2)=vp(:,2)/norm(vp(:,2));

H{1,1}=(1/sqrt(2))*[-0.9704 + 0.4012i (1.9144 - 0.3561i);(0.4516 - 1.4800i) 0.0501 - 0.1627i];
H{1,2}=(1/sqrt(2))*[-0.6337 + 0.8001i (-1.1485 + 0.2689i);(0.4516 - 1.4800i) -0.6651 - 0.9268i];
H{2,1}=(1/sqrt(2))*[-0.0165 + 0.9251i (-1.0463 - 0.5763i);(-0.1497 - 1.5829i) -0.7804 - 0.4109i];
H{2,2}=(1/sqrt(2))*[-0.9313 + 0.8060i (0.7313 + 0.0698i);(-0.2850 + 1.1345i) -0.5113 - 0.1662i];

sigma = sqrt(10^(-3));

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

            
         
     
                  
            gc_wiener(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}' + H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}' + sum_c1_f(:,k)*sum_c1_f(:,k)' +sum_p1_f(:,k)*sum_p1_f(:,k)' + H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2  ) * ( sum_c2_f(:,k) );
            gc_wiener(:,k) = gc_wiener(:,k)/norm(gc_wiener(:,k));
            gp_wiener(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}' + H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}' + sum_c1_f(:,k)*sum_c1_f(:,k)' +sum_p1_f(:,k)*sum_p1_f(:,k)' + H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2  ) * ( H{k,k}*vp(:,k) );
            gp_wiener(:,k) = gp_wiener(:,k)/norm(gp_wiener(:,k));
           
end

for k = 1:2 

            sum_c1_b(:,k) = [0;0];
            for j = 1:2
                if j~=k
                    sum_c1_b(:,k) = sum_c1_b(:,k) + H{j,k}.'*gc_wiener(:,j);
                end
            end

            sum_c2_b(:,k) = [0;0];
            for j = 1:2
                    sum_c2_b(:,k) = sum_c2_b(:,k) + H{j,k}.'*gc_wiener(:,j);
            end

            sum_p1_b(:,k) = [0;0];
            for j = 1:2
                if j~=k
                    sum_p1_b(:,k) = sum_p1_b(:,k) + H{j,k}.'*gp_wiener(:,j);
                end
            end
     
                  
            vc_wiener(:,k) = inv(  H{k,k}.'*gc_wiener(:,k)*gc_wiener(:,k)'*(H{k,k}').' + H{k,k}.'*gp_wiener(:,k)*gp_wiener(:,k)'*(H{k,k}').' + sum_c1_b(:,k)*sum_c1_b(:,k)' +sum_p1_b(:,k)*sum_p1_b(:,k)' + H{k,k}.'*gc_wiener(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc_wiener(:,k)'*(H{k,k}').'+eye(2)*sigma^2  ) * ( sum_c2_b(:,k) ); 
            vc_wiener(:,k) = vc_wiener(:,k)/norm(vc_wiener(:,k));
            vp_wiener(:,k) = inv(  H{k,k}.'*gc_wiener(:,k)*gc_wiener(:,k)'*(H{k,k}').' + H{k,k}.'*gp_wiener(:,k)*gp_wiener(:,k)'*(H{k,k}').' + sum_c1_b(:,k)*sum_c1_b(:,k)' +sum_p1_b(:,k)*sum_p1_b(:,k)' + H{k,k}.'*gc_wiener(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc_wiener(:,k)'*(H{k,k}').'+eye(2)*sigma^2  ) * ( H{k,k}.'*gp_wiener(:,k) );
            vp_wiener(:,k) = vp_wiener(:,k)/norm(vp_wiener(:,k));
            
            
end

SINR_C(1)=(   norm( gc_wiener(:,1)'*H{1,1}*vc_wiener(:,1)+gc_wiener(:,2)'*H{1,2}*vc_wiener(:,2) )^2   )/( norm( gc_wiener(:,1)'*sigma^2*gc_wiener(:,1) ) +  norm( gc_wiener(:,1)'*H{1,1}*vp_wiener(:,1)+gc_wiener(:,1)'*H{1,2}*vp_wiener(:,2) )^2  );
SINR_C(2)=(   norm( gc_wiener(:,2)'*H{2,1}*vc_wiener(:,1)+gc_wiener(:,2)'*H{2,2}*vc_wiener(:,2) )^2   )/( norm( gc_wiener(:,2)'*sigma^2*gc_wiener(:,2) ) +  norm( gc_wiener(:,2)'*H{2,1}*vp_wiener(:,1)+gc_wiener(:,2)'*H{2,2}*vp_wiener(:,2) )^2  );

SINR_P(1)=(   norm( gp_wiener(:,1)'*H{1,1}*vp_wiener(:,1) )^2   )/( norm( gp_wiener(:,1)'*sigma^2*gp_wiener(:,1) ) +  norm( gp_wiener(:,1)'*H{1,1}*vc_wiener(:,1)+gp_wiener(:,1)'*H{1,2}*vc_wiener(:,2)+gp_wiener(:,1)'*H{1,2}*vp_wiener(:,2) )^2  );
SINR_P(2)=(   norm( gp_wiener(:,2)'*H{2,2}*vp_wiener(:,2) )^2   )/( norm( gp_wiener(:,2)'*sigma^2*gp_wiener(:,2) ) +  norm( gp_wiener(:,2)'*H{2,1}*vc_wiener(:,1)+gp_wiener(:,2)'*H{2,2}*vc_wiener(:,2)+gp_wiener(:,2)'*H{2,1}*vp_wiener(:,1) )^2  );

C2=0;
for i = 1:2
    
    C2=C2+log2(1+SINR_C(i))+log2(1+SINR_P(i))
    
end