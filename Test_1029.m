%2 user, 2X2 MIMO Channel
%SINR(Stochastic Algorithm)
%calculate gc(:,1) gp(:,1) gc(:,2) gp(:,2) then V with cooperation

function [C] = Test_1029(times,sigma,StepSize,k,z)





for length = 1:z
    
    length
    
    C(length)=0;

    for zz = 1:times

        vc(:,1)=[1.4291 - 0.6656i;0.0183 - 1.1054i];
        vc(:,2)=[0.8221 - 0.3290i;0.3581 + 0.1167i];
        vp(:,1)=[1.9629 + 0.1063i;2.0188 + 0.3712i];
        vp(:,2)=[-0.0221 - 1.3322i;-0.4321 + 1.5764i];

        vc(:,1)=vc(:,1)/norm(vc(:,1));
        vc(:,2)=vc(:,2)/norm(vc(:,2));
        vp(:,1)=vp(:,1)/norm(vp(:,1));
        vp(:,2)=vp(:,2)/norm(vp(:,2));


        gc(:,1) = [-0.9072 + 0.2915i;-1.1637 - 0.6096i];
        gp(:,1) = [1.1707 - 0.0211i;-0.6295 - 1.1422i];
        gc(:,2) = [1.0489 - 0.1224i;-1.7782 - 2.3879i];
        gp(:,2) = [0.1764 - 1.7156i;0.8745 + 0.4395i];

        gc(:,1)=gc(:,1)/norm(gc(:,1));
        gc(:,2)=gc(:,2)/norm(gc(:,2));
        gp(:,1)=gp(:,1)/norm(gp(:,1));
        gp(:,2)=gp(:,2)/norm(gp(:,2));

        H{1,1}=(1/sqrt(2))*[-0.9704 + 0.4012i (1.9144 - 0.3561i);(0.4516 - 1.4800i) 0.0501 - 0.1627i];
        H{1,2}=0.8*(1/sqrt(2))*[-0.6337 + 0.8001i (-1.1485 + 0.2689i);(0.4516 - 1.4800i) -0.6651 - 0.9268i];
        H{2,1}=0.8*(1/sqrt(2))*[-0.0165 + 0.9251i (-1.0463 - 0.5763i);(-0.1497 - 1.5829i) -0.7804 - 0.4109i];
        H{2,2}=(1/sqrt(2))*[-0.9313 + 0.8060i (0.7313 + 0.0698i);(-0.2850 + 1.1345i) -0.5113 - 0.1662i];
        
        V = [vc(:,1);vc(:,2)];
        
        A{1,1}=H{1,1}.';
        A{1,2}=H{2,1}.';
        A{2,1}=H{1,2}.';
        A{2,2}=H{2,2}.';
        
        Z=[A{1,1} A{1,2};A{2,1} A{2,2}];

        Ah{1,1}=conj(H{1,1});
        Ah{1,2}=conj(H{1,2});
        Ah{2,1}=conj(H{2,1});
        Ah{2,2}=conj(H{2,2});
        
        Zh=[Ah{1,1} Ah{1,2};Ah{2,1} Ah{2,2}];

        x = zeros(1,length);
        xp = zeros(length,1);
        xp = zeros(length,2);

        for iteration = 1:k
        
            for iter = 1:length

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

                    gc(:,k) = gc(:,k)+StepSize*u(:,k)*conj(x(iter)-gc(:,k)'*u(:,k));
                    gp(:,k) = gp(:,k)+StepSize*u(:,k)*conj(xp(iter,k)-gp(:,k)'*u(:,k));

                end

            end

            for k = 1:2

            gc(:,k) = gc(:,k)/norm(gc(:,k));
            gp(:,k) = gp(:,k)/norm(gp(:,k));

            end

            for iter = 1:length

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


                    V_wiener = inv( Z*[gc(:,1);gc(:,2)]*[gc(:,1);gc(:,2)]'*Zh+Z*[gp(:,1)*gp(:,1)' zeros(2);zeros(2) gp(:,2)*gp(:,2)']*Zh +eye(4)*sigma^2    )*Z*[gc(:,1);gc(:,2)];
                    vc_wiener(:,k) = sqrt(2)*[V_wiener(k,1);V_wiener(k+1,1)]/norm(V_wiener);
                    vp_wiener(:,k) = inv(  H{k,k}.'*gc(:,k)*gc(:,k)'*(H{k,k}').' + H{k,k}.'*gp(:,k)*gp(:,k)'*(H{k,k}').' + sum_c1_b(:,k)*sum_c1_b(:,k)' +sum_p1_b(:,k)*sum_p1_b(:,k)' + H{k,k}.'*gc(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc(:,k)'*(H{k,k}').'+eye(2)*sigma^2  ) * ( H{k,k}.'*gp(:,k) );
                    vp_wiener(:,k) = vp_wiener(:,k)/norm(vp_wiener(:,k));

                    vp(:,k) = vp(:,k)+StepSize*u(:,k)*0.5*conj(xp(iter,k)-vp(:,k)'*u(:,k));
                    vp(:,k) = vp(:,k)+StepSize*u(:,k)*0.5*conj(xp(iter,k)-vp(:,k)'*u(:,k));


               end

               dummy(:,1) = vc(:,1);
               vc(:,1) = vc(:,1)+StepSize*u(:,1)*0.5*conj(x(iter)-vc(:,1)'*u(:,1)-vc(:,2)'*u(:,2));
               vc(:,2) = vc(:,2)+StepSize*u(:,2)*0.5*conj(x(iter)-vc(:,2)'*u(:,2)-dummy(:,1)'*u(:,1));
               dummy(:,1) = vc(:,1);
               vc(:,1) = sqrt(2)*vc(:,1)/norm([vc(:,1);vc(:,2)]);
               vc(:,2) = sqrt(2)*vc(:,2)/norm([dummy(:,1);vc(:,2)]);

            end

            for k = 1:2

                vp(:,k) = vp(:,k)/norm(vp(:,k));

            end
            
        end


        SINR_C(1)=(   norm( gc(:,1)'*H{1,1}*vc(:,1)+gc(:,2)'*H{1,2}*vc(:,2) )^2   )/( norm( gc(:,1)'*sigma^2*gc(:,1) ) +  norm( gc(:,1)'*H{1,1}*vp(:,1)+gc(:,1)'*H{1,2}*vp(:,2) )^2  );
        SINR_C(2)=(   norm( gc(:,2)'*H{2,1}*vc(:,1)+gc(:,2)'*H{2,2}*vc(:,2) )^2   )/( norm( gc(:,2)'*sigma^2*gc(:,2) ) +  norm( gc(:,2)'*H{2,1}*vp(:,1)+gc(:,2)'*H{2,2}*vp(:,2) )^2  );
        SINR_P(1)=(   norm( gp(:,1)'*H{1,1}*vp(:,1) )^2   )/( norm( gp(:,1)'*sigma^2*gp(:,1) ) +  norm( gp(:,1)'*H{1,1}*vc(:,1)+gp(:,1)'*H{1,2}*vc(:,2)+gp(:,1)'*H{1,2}*vp(:,2) )^2  );
        SINR_P(2)=(   norm( gp(:,2)'*H{2,2}*vp(:,2) )^2   )/( norm( gp(:,2)'*sigma^2*gp(:,2) ) +  norm( gp(:,2)'*H{2,1}*vc(:,1)+gp(:,2)'*H{2,2}*vc(:,2)+gp(:,2)'*H{2,1}*vp(:,1) )^2  );


        C(length)=C(length)+( log2(1+SINR_C(1))+log2(1+SINR_P(1))+log2(1+SINR_C(2))+log2(1+SINR_P(2)) )/times;

    end

end

plot(C)


end

