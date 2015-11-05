%1 user, 2X2 MIMO Channel
%Turn Off Private Channel 
%calculate MSE
%LS Filter
%Only Forward Direction
clc
clear

H=[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];

v=[1;1];
v=v/norm(v);
v_w=[1;1];
v_w=v_w/norm(v_w);
g=[1;1];
g_w=[1;1];
g=g/norm(g);
g_w=g_w/norm(g_w);

sigma = sqrt(10^(-3));
StepSize = 10^(-3);

R=H*H'+sigma^2*eye(2);
2/max(eig(R));


MSE = zeros(1,10^(7));
LS_LMS = zeros(1,10^(7));
MMSE = zeros(1,10^(7));
C = zeros(1,10^(7));
C_Wiener = zeros(1,10^(7));
SINR = zeros(1,10^(7));
SINR_C_Wiener= zeros(1,10^(7));



    
    for i=1:5
           

            %Forward Training  
            for iter2 = 1:i

                    if rand-0.5 >= 0
                                xf(iter2) = 1;
                            else
                                xf(iter2) = -1;
                    end

                    yf(:,iter2) = H*( v*xf(iter2) )+ sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))]; 
                    
                   
            end 
    
            g_w = inv(H*v_w*v_w'*H'+eye(2)*sigma^2)*H*v_w;
            g  = inv(yf*yf')*yf*xf';
            
            MSE(i) = real(  1-v'*H'*g-(v'*H'*g)'+g'*eye(2)*(sigma^2)*g+(v'*H'*g)'*(v'*H'*g) );
            SINR_without_stat(i)= norm(( g'*H*v ))^2/norm( g'*eye(2)*sigma^2*g ); 
            SINR_know_stat(i)= norm(( g_w'*H*v_w ))^2/norm( g_w'*eye(2)*sigma^2*g_w ); 
            MMSE(i) = real(   1-v'*H'* inv(H*v*v'*H'+eye(2)*sigma^2)*H*v);
                   
    end
    
   
          

   

n=1:i;

subplot(2,1,1)
plot(n,MSE(n))
legend('MSE')
xlabel('Iteration')
ylabel('MSE')
title('LS;1 User;Fixed 2X2 MIMO;Pilot Length=50;\mu=10^{-3}')
axis([1 i 0 10^(-2)])

subplot(2,1,2)
plot(   n,log2(1+SINR_without_stat(n)),n,log2(1+SINR_know_stat(n)))
legend('C(Bi-Directional Training)','C(Max-SINR)')
xlabel('Iteration')
ylabel('C')
title('LS;1 User;Fixed 2X2 MIMO;Pilot Length=50;\mu=10^{-3}')
axis([1 i 0 15])

