%1 user, 2X2 MIMO Channel
%Turn Off Private Channel 
%calculate MSE
%LS Filter
clc
clear

H=[-0.9704 + 0.4012i -0.4445 + 0.8804i;1.-0.7016 + 1.0288i -0.2290 - 0.3583i];

v=[1;1];
v=v/norm(v);

g=[1;1];
g=g/norm(g)


sigma = sqrt(10^(-3));
StepSize = 10^(-3);

R=H*H'+sigma^2+eye(2);
2/max(eig(R));


SINR = zeros(1,10^(7));
MSE = zeros(1,10^(7));
LS_LMS = zeros(1,10^(7));
MMSE = zeros(1,10^(7));
C = zeros(1,10^(7));
C_Wiener = zeros(1,10^(7));
SINR_C_Wiener= zeros(1,10^(7));


v_w = inv(H.'*g*g'*(H.')'+eye(2)*sigma^2)*H.'*g;
v_w = v_w/norm(v_w);
g_w = inv(H*v_w*v_w'*H'+eye(2)*sigma^2)*H*v_w;
SINR_w= norm( g_w'*H*v_w )^2/norm( g_w'*sigma^2*g_w ); 



for i = 1:100
    i
    
    REALIZATION = 10;
    for R = 1:REALIZATION

         v=1;
         v=v/norm(v);
         g=[1;1];
         g=g/norm(g);
      

        %Backward Training
        for iter1 = 1:i

                if rand-0.5 >= 0
                            xb(iter1) = 1;
                        else
                            xb(iter1) = -1;
                end

                yb(:,iter1) = H.'*g*xb(iter1)+sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)]; 
                
        end
   
        v  = inv(yb*yb')*yb*xb';
        v=v/norm(v);
        

        %Forward Training
    
        for iter2 = 1:i

                if rand-0.5 >= 0
                        xf(iter2) = 1;
                    else
                        xf(iter2) = -1;
                end

                yf(:,iter2) = H*( v*xf(iter2) )+ sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))]; 
                
        end 
        
        g = inv(yf*yf')*yf*xf';
        
        MSE(i) = real( MSE(i)+(1-v'*H'*g-(v'*H'*g)'+g'*g*(sigma^2)+(v'*H'*g)'*(v'*H'*g))/REALIZATION);
        SINR(i)= real( SINR(i)+(norm( g'*H*v/norm(v) )^2/norm( g'*sigma^2*g ))/REALIZATION); 
        MMSE(i) = real(   1-v'*H'* inv(H*v*v'*H'+eye(2)*sigma^2)*H*v);
        
        g = g/norm(g);
        
        
    end
    
end



n=1:i;

subplot(2,1,1)
plot(n,MSE(n))
legend('MSE(LMS)')
xlabel('Time n')
ylabel('MSE')
title('1 User;2X2 MIMO')

subplot(2,1,2)
plot(n,log2(1+SINR(n)),n,log2(1+SINR_w)+n-n)
legend('C')
xlabel('Time n')
ylabel('SINR')
title('1 User;2X2 MIMO')


