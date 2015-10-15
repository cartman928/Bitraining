clc
clear
%Beamformer for the very first iteration
vc(:,1) = [1;1]/norm([1;1]);
vp(:,1) = [1;1]/norm([1;1]);
vc(:,2) = [1;1]/norm([1;1]);
vp(:,2) = [1;1]/norm([1;1]);

%{
Channel
H{1,1}=[0.85 0.3;0.9 0.35];
H{1,2}=[0.7 0.7;0.55 0.75];
H{2,1}=[0.9 0.5;0.73 0.6];
H{2,2}=[0.69 0.73;0.69 0.91];
%}

H{1,1}=[1 1;1 1];
H{1,2}=[1 1;1 1];
H{2,1}=[1 1;1 1];
H{2,2}=[1 1;1 1];

sigma = 10^(-2);

for iter = 1:1000 %Number of iterations

        %Random Source Message
        %Common Message
        if rand-0.5 >= 0
                    x(iter) = 1;
                else
                    x(iter) = -1;
        end

        %Private Message for User#1
        if rand-0.5 >= 0
                    xp_1(iter) = 1;
                else
                    xp_1(iter) = -1;
        end

        %Private Message for User#2
        if rand-0.5 >= 0
                    xp_2(iter) = 1;
                else
                    xp_2(iter) = -1;
        end


        %Calculate Receivers then Transmitters

        for k = 1:2 %For user 1 and 2

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


            gc(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}' + H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}' + sum_c1_f(:,k)*sum_c1_f(:,k)' +sum_p1_f(:,k)*sum_p1_f(:,k)' + H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2  ) * ( sum_c2_f(:,k) );     
            gc(:,k) = gc(:,k)/norm(gc(:,k))%normailize
            gp(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}' + H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}' + sum_c1_f(:,k)*sum_c1_f(:,k)' +sum_p1_f(:,k)*sum_p1_f(:,k)' + H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2  ) * ( H{k,k}*vp(:,k) );
            gp(:,k) = gp(:,k)/norm(gp(:,k))%normailize

        end
        
        %%Decoding
        bc_1(iter)=sign(  gc(:,1)'*(      H{1,1}*(x(iter)*vc(:,1)+xp_1(iter)*vp(:,1))+H{1,2}*(x(iter)*vc(:,2)+xp_2(iter)*vp(:,2)) + sigma*[rand;rand] )  );
        bp_1(iter)=sign(  gp(:,1)'*(      H{1,1}*(x(iter)*vc(:,1)+xp_1(iter)*vp(:,1))+H{1,2}*(x(iter)*vc(:,2)+xp_2(iter)*vp(:,2)) + sigma*[rand;rand] )  );
        bc_2(iter)=sign(  gc(:,2)'*(      H{2,1}*(x(iter)*vc(:,1)+xp_1(iter)*vp(:,1))+H{2,2}*(x(iter)*vc(:,2)+xp_2(iter)*vp(:,2)) + sigma*[rand;rand] )  );
        bp_2(iter)=sign(  gp(:,2)'*(      H{2,1}*(x(iter)*vc(:,1)+xp_1(iter)*vp(:,1))+H{2,2}*(x(iter)*vc(:,2)+xp_2(iter)*vp(:,2)) + sigma*[rand;rand] )  );


        for k = 1:2 %For user 1 and 2

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

            vc(:,k) = inv(  H{k,k}.'*gc(:,k)*gc(:,k)'*(H{k,k}.')'+sum_c1_b(:,k)*sum_c1_b(:,k)'+H{k,k}.'*gc(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc(:,k)'*(H{k,k}.')'+eye(2)*sigma^2+H{k,k}.'*gp(:,k)*gp(:,k)'*(H{k,k}.')'+sum_p1_b(:,k)*sum_p1_b(:,k)'  ) * ( sum_c2_b(:,k) );
            vc(:,k) = vc(:,k)/norm(vc(:,k));%normailize
            vp(:,k) = inv(  H{k,k}.'*gc(:,k)*gc(:,k)'*(H{k,k}.')'+sum_c1_b(:,k)*sum_c1_b(:,k)'+H{k,k}.'*gc(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc(:,k)'*(H{k,k}.')'+eye(2)*sigma^2+H{k,k}.'*gp(:,k)*gp(:,k)'*(H{k,k}.')'+sum_p1_b(:,k)*sum_p1_b(:,k)'  ) * ( H{k,k}.'*gp(:,k) );
            vp(:,k) = vp(:,k)/norm(vp(:,k));%normailize
        end

end

%{
error=0;
for k = 101:1000
                    error = error + abs(x(k)-bc_1(k))+abs(x(k)-bc_2(k))+abs(x(k)-bc_2(k))+abs(xp_1(k)-bp_1(k))+abs(xp_2(k)-bp_2(k));
            end
4*900-error
%}

%{
x(9991:10000)
bc_1(9991:10000)
bc_2(9991:10000)

xp_1(9991:10000)
bp_1(9991:10000)
xp_2(9991:10000)
bp_2(9991:10000)

x(91:100)
bc_1(91:100)
bc_2(91:100)

xp_1(91:100)
bp_1(91:100)
xp_2(91:100)
bp_2(91:100)
%}


 
