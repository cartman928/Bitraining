%Beam Steering
%Choose Starting Point to Be Zero Vector!!!!!!

clc
clear

H=[-0.9704 + 0.4012i;1.0291 - 0.4917i];

w=[0;0];

sigma = sqrt(10^(-3));
StepSize = 10^(-3);

R=H*H'+sigma^2+eye(2);
e = eig(R);

w_wiener = inv(H*H'+sigma^2)*H;
SINR_wiener = norm(w_wiener'*H)^2/norm(w_wiener'*sigma^2*w_wiener); 

for iteration = 1:10^(2)
    
            iteration

            if rand-0.5 >= 0
                        x(iteration) = 1;
                    else
                        x(iteration) = -1;
            end

            y= H*x(iteration)+sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))]; 
            w = w+StepSize*y*conj(x(iteration)-w'*y);
            
            MSE(iteration) = real( 1-H'*w-w'*H+w'*H*H'*w+w'*(sigma^2)*eye(2)*w  );
            SNR(iteration)= norm( w'*H )^2/norm( w'*sigma^2*eye(2)*w ); 
            MMSE(iteration)= real( 1-H'*inv(H*H'+eye(2)*sigma^2)*H );
            
end 
    

   
    
x=1:iteration;

subplot(2,1,1)
plot(x,MSE(x),x,MMSE(x))
legend('MSE','MMSE')
xlabel('Time n')
ylabel('MSE')
title('Beam Steering;2 Antenna')

subplot(2,1,2)
plot(x,log2(1+SNR(x)))
legend('SINR')
xlabel('Time n')
ylabel('SINR')
title('Beam Steering;2 Antenna')


