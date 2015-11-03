%1 user, 2X2 MIMO Channel
%calculate G 
%calculate MSE
%Only 1 iteration
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

i = 10^(3); %FilterLength  
    
    %Backward Training
 
    for iter1 = 1:i

            if rand-0.5 >= 0
                        xb(iter1) = 1;
                    else
                        xb(iter1) = -1;
            end

            yb = H.'*g*xf(iter1)+sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)]; 
            v  = v+StepSize*yb*conj(xb(iter1)-v'*yb);

    end 
    
    %Normalize Transmitters
    v=v/norm(v);
    
    %Forward Training
    g=[0;0];
    for iter2 = 1:i
        
            iter2

            if rand-0.5 >= 0
                        xf(iter2) = 1;
                    else
                        xf(iter2) = -1;
            end

            yf = H*( v*xf(iter2) )+ sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))]; 
            g = g+StepSize*yf*conj(xf(iter2)-g'*yf);
            
            MSE(iter2) = 1-v'*H'*g-(v'*H'*g)'+g'*g*(sigma^2)+(v'*H'*g)'*(v'*H'*g);
            SINR(iter2)= norm( g'*H*v )^2/norm( g'*sigma^2*g ); 
            MMSE(iter2) = real(   1-v'*H'* inv(H*v*v'*H'+eye(2)*sigma^2)*H*v);
                   
    end 
    
    

n=1:i;

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


