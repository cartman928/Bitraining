clc
clear

%{
gc(:,1) = [1;1];
gp(:,1) = [1;1];
vc(:,1) = [1,1];
vp(:,1) = [1,1];
H{1,1}=(1/sqrt(2))*[-0.9704 + 0.4012i 0.8*(1.9144 - 0.3561i);0.8*(0.4516 - 1.4800i) 0.0501 - 0.1627i];
%}


gc(:,1) = 1;
gp(:,1) = 1;
vc(:,1) = 1;
vp(:,1) = 1;
H{1,1}=-0.9704 + 0.4012i;


sigma = 10^(-7);
StepSize = 10^(-4)

for iter = 1:10^(5) 
        if rand-0.5 >= 0
                    x(iter) = 1;
                else
                    x(iter) = -1;
        end

        if rand-0.5 >= 0
                    xp_1(iter) = 1;
                else
                    xp_1(iter) = -1;
        end


gc_wiener(:,1) = inv(  H{1,1}*vc(:,1)*vc(:,1)'*H{1,1}' + H{1,1}*vp(:,1)*vp(:,1)'*H{1,1}' +eye(1)*sigma^2  ) * ( H{1,1}*vc(:,1) );    
gp_wiener(:,1) = inv(  H{1,1}*vc(:,1)*vc(:,1)'*H{1,1}' + H{1,1}*vp(:,1)*vp(:,1)'*H{1,1}' +eye(1)*sigma^2  ) * ( H{1,1}*vp(:,1) );  

u(:,1) =  H{1,1}*( vc(:,1)*x(iter)+vp(:,1)*xp_1(iter) ) + sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1)];

gc(:,1) = gc(:,1)+StepSize*u(:,1)*conj(x(iter)-gc(:,1)'*u(:,1))
gp(:,1) = gp(:,1)+StepSize*u(:,1)*conj(xp_1(iter)-gp(:,1)'*u(:,1))

end

gc_wiener(:,1)
gp_wiener(:,1)
