%1 user, 4X4 MIMO Channel
%calculate G 
%calculate MSE
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


G = [gc(:,1);gc(:,2)];

H{1,1}=(1/sqrt(2))*[-0.9704 + 0.4012i (1.9144 - 0.3561i);(0.4516 - 1.4800i) 0.0501 - 0.1627i];
H{1,2}=0.8*(1/sqrt(2))*[-0.6337 + 0.8001i (-1.1485 + 0.2689i);(0.4516 - 1.4800i) -0.6651 - 0.9268i];
H{2,1}=0.8*(1/sqrt(2))*[-0.0165 + 0.9251i (-1.0463 - 0.5763i);(-0.1497 - 1.5829i) -0.7804 - 0.4109i];
H{2,2}=(1/sqrt(2))*[-0.9313 + 0.8060i (0.7313 + 0.0698i);(-0.2850 + 1.1345i) -0.5113 - 0.1662i];

B=[H{1,1} H{1,2};H{2,1} H{2,2}];


sigma = sqrt(10^(-3));
StepSize = 10^(-6);
x = zeros(1,10^(7));
xp(:,1) = zeros(1,10^(7));
xp(:,2) = zeros(1,10^(7));
MSE_LMS = zeros(1,10^(7));
MSE_Wiener = zeros(1,10^(7));
C_LMS = zeros(1,10^(7));
C_Wiener = zeros(1,10^(7));

G_Wiener = inv( B*[vc(:,1);vc(:,2)]*[vc(:,1);vc(:,2)]'*B'+ eye(4)*sigma^2    )*B*[vc(:,1);vc(:,2)];

for iter = 1:10^(5) 
    
    iter
    
        if rand-0.5 >= 0
                    x(iter) = 1;
                else
                    x(iter) = -1;
        end

        

            u(:,1) = [0;0];
            u(:,2) = [0;0];
            for j = 1:2
                    u(:,1) = u(:,1) + H{1,j}*( vc(:,j)*x(iter) );
                    u(:,2) = u(:,2) + H{2,j}*( vc(:,j)*x(iter) );
            end
            u(:,1) = u(:,1)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
            u(:,2) = u(:,2)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];

    
           
       G = G+StepSize*[u(:,1);u(:,2)]*conj(x(iter)-G'*[u(:,1);u(:,2)]);
       
       MSE_LMS(iter) = (G'*[u(:,1);u(:,2)]- x(iter))'*(G'*[u(:,1);u(:,2)]- x(iter));
       MSE_Wiener(iter) = (G_Wiener'*[u(:,1);u(:,2)]- x(iter))'*(G_Wiener'*[u(:,1);u(:,2)]- x(iter));
       
      
       SINR_C(1)=(   norm( [G(1,1);G(2,1)]'*H{1,1}*vc(:,1)+[G(1,1);G(2,1)]'*H{1,2}*vc(:,2) )^2   )/( norm( [G(1,1);G(2,1)]'*sigma^2*[G(1,1);G(2,1)] )    );
       SINR_C(2)=(   norm( [G(3,1);G(4,1)]'*H{2,1}*vc(:,1)+[G(3,1);G(4,1)]'*H{2,2}*vc(:,2) )^2   )/( norm( [G(3,1);G(4,1)]'*sigma^2*[G(3,1);G(4,1)] )    );
       C_LMS(iter)=log2(1+SINR_C(1))+log2(1+SINR_C(2)); 
   
       SINR_C_Wiener(1)=(   norm( [G_Wiener(1,1);G_Wiener(2,1)]'*H{1,1}*vc(:,1)+[G_Wiener(1,1);G_Wiener(2,1)]'*H{1,2}*vc(:,2) )^2   )/( norm( [G_Wiener(1,1);G_Wiener(2,1)]'*sigma^2*[G_Wiener(1,1);G_Wiener(2,1)] )    );
       SINR_C_Wiener(2)=(   norm( [G_Wiener(3,1);G_Wiener(4,1)]'*H{2,1}*vc(:,1)+[G_Wiener(3,1);G_Wiener(4,1)]'*H{2,2}*vc(:,2) )^2   )/( norm( [G_Wiener(3,1);G_Wiener(4,1)]'*sigma^2*[G_Wiener(3,1);G_Wiener(4,1)] )    );
       C_Wiener(iter)=log2(1+SINR_C_Wiener(1))+log2(1+SINR_C_Wiener(2));
       


end



x=1:iter;


subplot(2,1,1)
plot(x,MSE_LMS(x),x,MSE_Wiener(x))
legend('MSE(LMS)','MSE(Wiener)')
xlabel('Time n')
ylabel('MSE')
title('1User;4X4 MIMO;Only Common Msg;M=1')

subplot(2,1,2)
plot(x,C_LMS(x),x,C_Wiener(x))
legend('C(LMS)','C(Wiener)')
xlabel('Time n')
ylabel('C(bit/channel)')
title('1User;4X4 MIMO;Only Common Msg;M=1')



