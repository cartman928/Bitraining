%2 user, 2X2 MIMO Channel
%SINR
clc
clear

gc(:,1) = [0;0];
gp(:,1) = [0;0];
gc(:,2) = [0;0];
gp(:,2) = [0;0];

vp(:,1) = [1.0562 + 0.2200i,0.4991 - 0.4516i];
vp(:,1) = vp(:,1)/norm(vp(:,1))
vp(:,2) = [-1.3609 + 0.8642i,0.2347 - 0.0695i];
vp(:,2) = vp(:,2)/norm(vp(:,2));

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

            
         
     
                  
            %gc_wiener(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}' + H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}' + sum_c1_f(:,k)*sum_c1_f(:,k)' +sum_p1_f(:,k)*sum_p1_f(:,k)' + H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2  ) * ( sum_c2_f(:,k) );
            %gc_wiener(:,k) = gc_wiener(:,k)/norm(gc_wiener(:,k));
            gp_wiener(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}' + H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}' + sum_c1_f(:,k)*sum_c1_f(:,k)' +sum_p1_f(:,k)*sum_p1_f(:,k)' + H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2  ) * ( H{k,k}*vp(:,k) );
            gp_wiener(:,k) = gp_wiener(:,k)/norm(gp_wiener(:,k));
            gc_wiener(:,k) = [0;0];
            
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
     
                  
            %vc_wiener(:,k) = inv(  H{k,k}.'*gc_wiener(:,k)*gc_wiener(:,k)'*(H{k,k}').' + H{k,k}.'*gp_wiener(:,k)*gp_wiener(:,k)'*(H{k,k}').' + sum_c1_b(:,k)*sum_c1_b(:,k)' +sum_p1_b(:,k)*sum_p1_b(:,k)' + H{k,k}.'*gc_wiener(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc_wiener(:,k)'*(H{k,k}').'+eye(2)*sigma^2  ) * ( sum_c2_b(:,k) ); 
            %vc_wiener(:,k) = vc_wiener(:,k)/norm(vc_wiener(:,k));
            vp_wiener(:,k) = inv(  H{k,k}.'*gc_wiener(:,k)*gc_wiener(:,k)'*(H{k,k}').' + H{k,k}.'*gp_wiener(:,k)*gp_wiener(:,k)'*(H{k,k}').' + sum_c1_b(:,k)*sum_c1_b(:,k)' +sum_p1_b(:,k)*sum_p1_b(:,k)' + H{k,k}.'*gc_wiener(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc_wiener(:,k)'*(H{k,k}').'+eye(2)*sigma^2  ) * ( H{k,k}.'*gp_wiener(:,k) );
            vp_wiener(:,k) = vp_wiener(:,k)/norm(vp_wiener(:,k));
            vc_wiener(:,k) = [0;0];
            
end

%SINR_C(1)=(   norm( gc_wiener(:,1)'*H{1,1}*vc_wiener(:,1)+gc_wiener(:,2)'*H{1,2}*vc_wiener(:,2) )^2   )/( norm( gc_wiener(:,1)'*sigma^2*gc_wiener(:,1) ) +  norm( gc_wiener(:,1)'*H{1,1}*vp_wiener(:,1)+gc_wiener(:,1)'*H{1,2}*vp_wiener(:,2) )^2  );
%SINR_C(2)=(   norm( gc_wiener(:,2)'*H{2,1}*vc_wiener(:,1)+gc_wiener(:,2)'*H{2,2}*vc_wiener(:,2) )^2   )/( norm( gc_wiener(:,2)'*sigma^2*gc_wiener(:,2) ) +  norm( gc_wiener(:,2)'*H{2,1}*vp_wiener(:,1)+gc_wiener(:,2)'*H{2,2}*vp_wiener(:,2) )^2  );

SINR_P(1)=(   norm( gp_wiener(:,1)'*H{1,1}*vp_wiener(:,1) )^2   )/( norm( gp_wiener(:,1)'*sigma^2*gp_wiener(:,1) ) +  norm( gp_wiener(:,1)'*H{1,1}*vc_wiener(:,1)+gp_wiener(:,1)'*H{1,2}*vc_wiener(:,2)+gp_wiener(:,1)'*H{1,2}*vp_wiener(:,2) )^2  );
SINR_P(2)=(   norm( gp_wiener(:,2)'*H{2,2}*vp_wiener(:,2) )^2   )/( norm( gp_wiener(:,2)'*sigma^2*gp_wiener(:,2) ) +  norm( gp_wiener(:,2)'*H{2,1}*vc_wiener(:,1)+gp_wiener(:,2)'*H{2,2}*vc_wiener(:,2)+gp_wiener(:,2)'*H{2,1}*vp_wiener(:,1) )^2  );

C2=0;
for i = 1:2
    
    C2=C2+log2(1+SINR_P(i))+log2(1+SINR_C(i))
    
end