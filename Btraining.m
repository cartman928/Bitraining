%Beamformer for the very first iteration

vc(:,1) = [1;1];
vp(:,1) = [1;1];
vc(:,2) = [1;1];
vp(:,2) = [1;1];

%Channel

H{1,1}=[0.8 0.7;0.9 0.75];
H{1,2}=[0.7 0.7;0.75 0.75];
H{2,1}=[0.9 0.93;0.73 0.6];
H{2,2}=[0.69 0.73;0.69 0.91];
sigma = 0.001;%Square Root of Noise Variance

for i = 1:10 %Number of iterations

        %Random Source Message
        %Common Message
        if rand-0.5 >= 0
                    x(i) = 1;
                else
                    x(i) = -1;
        end

        %Private Message for User#1
        if rand-0.5 >= 0
                    xp_1(i) = 1;
                else
                    xp_1(i) = -1;
        end

        %Private Message for User#2
        if rand-0.5 >= 0
                    xp_2(i) = 1;
                else
                    xp_2(i) = -1;
        end


        %Calculate Receivers then Transmitters

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


            gc(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}'+sum_c1_f(:,k)*sum_c1_f(:,k)'+H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2+H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}'+sum_p1_f(:,k)*sum_p1_f(:,k)'  ) * ( sum_c2_f(:,k) );
            gp(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}'+sum_c1_f(:,k)*sum_c1_f(:,k)'+H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2+H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}'+sum_p1_f(:,k)*sum_p1_f(:,k)'  ) * ( H{k,k}*vp(:,k) )

        end


        for k = 1:2

            sum_c1_b(:,k) = [0;0];
            for j = 1:2
                if j~=k
                    sum_c1_b(:,k) = sum_c1_b(:,k) + H{k,j}.'*gc(:,j);
                end
            end

            sum_c2_b(:,k) = [0;0];
            for j = 1:2
                    sum_c2_b(:,k) = sum_c2_b(:,k) + H{k,j}.'*gc(:,j);
            end

            sum_p1_b(:,k) = [0;0];
            for j = 1:2
                if j~=k
                    sum_p1_b(:,k) = sum_p1_b(:,k) + H{k,j}.'*gp(:,j);
                end
            end

            vc(:,k) = inv(  H{k,k}.'*gc(:,k)*gc(:,k)'*(H{k,k}').'+sum_c1_b(:,k)*sum_c1_b(:,k)'+H{k,k}.'*gc(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc(:,k)'*(H{k,k}').'+eye(2)*sigma^2+H{k,k}.'*gp(:,k)*gp(:,k)'*(H{k,k}').'+sum_p1_b(:,k)*sum_p1_b(:,k)'  ) * ( sum_c2_b(:,k) );
            vp(:,k) = inv(  H{k,k}.'*gc(:,k)*gc(:,k)'*(H{k,k}').'+sum_c1_b(:,k)*sum_c1_b(:,k)'+H{k,k}.'*gc(:,k)*sum_c1_b(:,k)'+sum_c1_b(:,k)*gc(:,k)'*(H{k,k}').'+eye(2)*sigma^2+H{k,k}.'*gp(:,k)*gp(:,k)'*(H{k,k}').'+sum_p1_b(:,k)*sum_p1_b(:,k)'  ) * ( H{k,k}.'*gp(:,k) );

        end

end







%k=2;
%gc(k)=2;


%{
Decoded Messages

x1 = gc_1*RX1; 
xp_1=gp_1*RX1; 
x2 = gc_2*RX2; 
xp_2=gp_2*RX2;

%}



