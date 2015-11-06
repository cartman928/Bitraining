%2 user, 2X2 MIMO Channel
%Turn Off Private Channel 
%calculate MSE
clc
clear

H=[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];

v=[0;0];
g=[1;1];
g_w=[1;1]
g=g/norm(g);
g_w=g_w/norm(g_w);

sigma = sqrt(10^(-3));
StepSize = 10^(-3);

R=H*H'+sigma^2*eye(2);
2/max(eig(R));


xf = zeros(1,10^(7));
MSE = zeros(1,10^(7));
LS_LMS = zeros(1,10^(7));
MMSE = zeros(1,10^(7));
C = zeros(1,10^(7));
C_Wiener = zeros(1,10^(7));
SINR_know_stat = zeros(1,10^(7));
SINR_without_stat= zeros(1,10^(7));

Realization = 1;
TrainingLength = 20; %TrainingLength


    for R = 1:Realization
    
    R    
        
    v=[0;0];
    g=[1;1];
    g_w=[1;1];
    g=g/norm(g);
    g_w=g_w/norm(g_w);
    H=[randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1) randn(1,1)+1i*randn(1,1)];

        for iteration = 1:30

            iteration;

            for loop=1:iteration

                
                    g=g/norm(g);
                    g_w=g_w/norm(g_w);

                    %Backward Training
                    for iter1 = 1:TrainingLength

                            if rand-0.5 >= 0
                                        xb(iter1) = 1;
                                    else
                                        xb(iter1) = -1;
                            end

                            yb = H.'*g*xb(iter1)+sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)]; 
                       
                            v  = v+StepSize*yb*conj(xb(iter1)-v'*yb);
                            v_w = inv(H.'*g_w*g_w'*(H.')'+eye(2)*sigma^2)*H.'*g_w;

                    end

                    v=v/norm(v);
                    v_w=v_w/norm(v_w);

                    %Forward Training  
                    for iter2 = 1:TrainingLength

                            if rand-0.5 >= 0
                                        xf(iter2) = 1;
                                    else
                                        xf(iter2) = -1;
                            end

                            yf = H*( v*xf(iter2) )+ sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))]; 
                            g = g+StepSize*yf*conj(xf(iter2)-g'*yf);
                            g_w = inv(H*v_w*v_w'*H'+eye(2)*sigma^2)*H*v_w;
                      
                            
                            

            end 
            end

            MSE(iteration) = MSE(iteration) + real(  1-v'*H'*g-(v'*H'*g)'+g'*eye(2)*(sigma^2)*g+(v'*H'*g)'*(v'*H'*g) )/Realization;
            SINR_without_stat(iteration)= SINR_without_stat(iteration) + norm(( g'*H*v ))^2/norm( g'*eye(2)*sigma^2*g )/Realization; 
            SINR_know_stat(iteration)= SINR_know_stat(iteration)+norm(( g_w'*H*v_w ))^2/norm( g_w'*eye(2)*sigma^2*g_w )/Realization; 
            MMSE(iteration) = real(   1-v'*H'* inv(H*v*v'*H'+eye(2)*sigma^2)*H*v);
            
        end    
end
   

n=1:iteration;

subplot(2,1,1)
plot(n,MSE(n))
legend('MSE')
xlabel('Iteration')
ylabel('MSE')
title('1 User;2X2 MIMO;LMS')
axis([1 iteration 0 3])

subplot(2,1,2)
plot(   n,log2(1+SINR_without_stat(n)),n,log2(1+SINR_know_stat(n)))
legend('C(Bi-Directional Training)','C(Max-SINR)')
xlabel('Iteration')
ylabel('C')
title('1 User;2X2 MIMO;LMS')
axis([1 iteration 0 15])
