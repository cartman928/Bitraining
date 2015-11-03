%2 user, 2X2 MIMO Channel
%Turn Off Private Channel 
%calculate MSE
%Know Statistics
clc
clear

H=[-0.9704 + 0.4012i 0.2969 + 0.2337i;1.-0.7016 + 1.0288i 2.0200 - 0.1294i];

v=[0;0];

g=[1;1];
g=g/norm(g);

sigma = sqrt(10^(-3));
StepSize = 10^(-3);

R=H*H'+sigma^2+eye(2);
2/max(eig(R));


xf = zeros(1,10^(7));
MSE = zeros(1,10^(7));
LS_LMS = zeros(1,10^(7));
MMSE = zeros(1,10^(7));
C = zeros(1,10^(7));
C_Wiener = zeros(1,10^(7));
SINR = zeros(1,10^(7));
SINR_C_Wiener= zeros(1,10^(7));


%g_w = inv(H*v*v'*H'+sigma^2)*H*v;

i = 500; %FilterLength

for iteration = 1:1000

    iteration
    
    
    %Backward Training
    for iter1 = 1:i

            v  = inv(H.'*g*g'*(H.')'+sigma^2)*H.'*g;

    end
   
    v=v/norm(v);
    


    %Forward Training
    for iter2 = 1:i
  
            g = inv(H*v*v'*H'+sigma^2)*H*v;
                   
    end 
    
 %}
    
    MSE(iteration) = 1-v'*H'*g-(v'*H'*g)'+g'*g*(sigma^2)+(v'*H'*g)'*(v'*H'*g);
    SINR(iteration)= norm( g'*H*v/norm(v) )^2/norm( g'*sigma^2*g ); 
    MMSE(iteration) = real(   1-v'*H'* inv(H*v*v'*H'+eye(2)*sigma^2)*H*v);
 
    g=g/norm(g);   
        
end
   

n=1:iteration;

subplot(2,1,1)
plot(n,MSE(n))
legend('MSE(LMS)')
xlabel('Time n')
ylabel('MSE')
title('1 User;2X2 MIMO')

subplot(2,1,2)
plot(n,log2(1+SINR(n)))
legend('C')
xlabel('Time n')
ylabel('SINR')
title('1 User;2X2 MIMO')


