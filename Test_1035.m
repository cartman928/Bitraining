%1 user, 1X2 MIMO Channel
%Turn Off Common Channel 
%calculate MSE
%Statistics are known at receivers
clc
clear

H=[-0.9704 + 0.4012i ;1.-0.7016 + 1.0288i];

v=1;
v=v/norm(v);

g=[0;0];
%g=g/norm(g);

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

for iteration = 1:20

    iteration
    
    %{
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
    %}
    


    %Forward Training
    
    %if iteration == 1
    %g=[0;0];
    %end
    
    %for iter2 = 1:i

            if rand-0.5 >= 0
                        xf(iteration) = 1;
                    else
                        xf(iteration) = -1;
            end

            yf = H*( v*xf(iteration) )+ sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))]; 
            g = (eye(2)-StepSize*(H*H'+sigma^2*eye(2)))*g+StepSize*H;
                   
    %end 
    
 %}
    
    MSE(iteration) = 1-v'*H'*g-(v'*H'*g)'+g'*g*(sigma^2)+(v'*H'*g)'*(v'*H'*g);
    SINR(iteration)= norm( g'*H*v/norm(v) )^2/norm( g'*sigma^2*g ); 
    MMSE(iteration) = real(   1-v'*H'* inv(H*v*v'*H'+eye(2)*sigma^2)*H*v);
    
    %SINR_w= norm( g_w'*H*v )^2/norm( g_w'*sigma^2*g_w ); 
    %MMSE_w = real(  1-v'*H'* inv(H*v*v'*H'+eye(2)*sigma^2)*H*v);
    %g=g/norm(g);   
        
end
   

n=1:iteration;

subplot(2,1,1)
plot(n,MSE(n))
legend('MSE(LMS)')
xlabel('Time n')
ylabel('MSE')
title('1 User;1X2 MIMO')

subplot(2,1,2)
plot(n,log2(1+SINR(n)))
legend('C')
xlabel('Time n')
ylabel('SINR')
title('1 User;1X2 MIMO')


