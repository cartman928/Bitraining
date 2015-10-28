%2 user, 2X2 MIMO Channel
%calculate V,vp(:,1),vp(:,2) with cooperation
clc
clear


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

gp(:,1)=[0;0];
gp(:,2)=[0;0];

V = [vc(:,1);vc(:,2)];

H{1,1}=(1/sqrt(2))*[-0.9704 + 0.4012i (1.9144 - 0.3561i);(0.4516 - 1.4800i) 0.0501 - 0.1627i];
H{1,2}=0.8*(1/sqrt(2))*[-0.6337 + 0.8001i (-1.1485 + 0.2689i);(0.4516 - 1.4800i) -0.6651 - 0.9268i];
H{2,1}=0.8*(1/sqrt(2))*[-0.0165 + 0.9251i (-1.0463 - 0.5763i);(-0.1497 - 1.5829i) -0.7804 - 0.4109i];
H{2,2}=(1/sqrt(2))*[-0.9313 + 0.8060i (0.7313 + 0.0698i);(-0.2850 + 1.1345i) -0.5113 - 0.1662i];

B=[H{1,1} H{1,2};H{2,1} H{2,2}];

%{
Z{1,1}=H{1,1}.';
Z{1,2}=H{2,1}.';
Z{2,1}=H{1,2}.';
Z{2,2}=H{2,2}.';

Zh{1,1}=conj(H{1,1});
Zh{1,2}=conj(H{1,2});
Zh{2,1}=conj(H{2,1});
Zh{2,2}=conj(H{2,2});
%}


sigma = sqrt(10^(-3));
StepSize = 10^(-6);
x = zeros(1,10^(6));
xp(:,1) = zeros(1,10^(6));
xp(:,2) = zeros(1,10^(6));
MSE_LMS = zeros(1,10^(6));
MSE_Wiener = zeros(1,10^(6));

V_wiener = inv( B*[gc(:,1);gc(:,2)]*[gc(:,1);gc(:,2)]'*B'+B*[gp(:,1)*gp(:,1)' zeros(2);zeros(2) gp(:,2)*gp(:,2)']*B' +eye(4)*sigma^2    )*B*[gc(:,1);gc(:,2)];


for iter = 1:10^(6) 
    
    iter
    
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

            u(:,1) = [0;0];
            u(:,2) = [0;0];
            for j = 1:2
                    u(:,1) = u(:,1) + H{1,j}*( gc(:,j)*x(iter)+gp(:,j)*xp(iter,j) );
                    u(:,2) = u(:,2) + H{2,j}*( gc(:,j)*x(iter)+gp(:,j)*xp(iter,j) );
            end
            u(:,1) = u(:,1)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
            u(:,2) = u(:,2)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];

            %vp_wiener(:,k) = inv(  H{k,k}*vc(:,k)*vc(:,k)'*H{k,k}' + H{k,k}*vp(:,k)*vp(:,k)'*H{k,k}' + sum_c1_f(:,k)*sum_c1_f(:,k)' +sum_p1_f(:,k)*sum_p1_f(:,k)' + H{k,k}*vc(:,k)*sum_c1_f(:,k)'+sum_c1_f(:,k)*vc(:,k)'*H{k,k}'+eye(2)*sigma^2  ) * ( H{k,k}*vp(:,k) );
       
            %vp(:,k) = vp(:,k)+StepSize*u(:,k)*conj(xp(iter,k)-vp(:,k)'*u(:,k));

    
       
       
       V = V+StepSize*[u(:,1);u(:,2)]*conj(x(iter)-V'*[u(:,1);u(:,2)]);
       
       MSE_LMS(iter) = (V'*[u(:,1);u(:,2)]- x(iter))'*(V'*[u(:,1);u(:,2)]- x(iter));
       MSE_Wiener(iter) = (V_wiener'*[u(:,1);u(:,2)]- x(iter))'*(V_wiener'*[u(:,1);u(:,2)]- x(iter));
       

end


x=1:iter;
plot(x,MSE_LMS(x),x,MSE_Wiener(x))
legend('MSE_LMS','MSE_Wiener')


