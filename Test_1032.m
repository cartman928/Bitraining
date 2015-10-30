%1 user, 1X1 MIMO Channel
%calculate G 
%calculate MSE
clc
clear

H=[-0.9704 + 0.4012i;1.0291 - 0.4917i];

vc(:,1)=[-0.8051+1.3128i];
vc(:,1)=vc(:,1)/norm(vc(:,1));

gc(:,1)=[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))];
gc(:,1)=gc(:,1)/norm(gc(:,1));



%2/real(H*vc(:,1)*vc(:,1)'*H'+sigma^2)

sigma = sqrt(10^(-3));
StepSize = 10^(-3);

R=H*H'+sigma^2+eye(2);
e = eig(R)

%{
x = zeros(1,10^(7));
MSE_LMS = zeros(1,10^(7));
LS_LMS = zeros(1,10^(7));
MMSE = zeros(1,10^(7));
C = zeros(1,10^(7));
C_Wiener = zeros(1,10^(7));
SINR = zeros(1,10^(7));
SINR_C_Wiener= zeros(1,10^(7));
%}

gc_w(:,1) = inv( H*vc(:,1)*vc(:,1)'*H'+ sigma^2    )*H*vc(:,1);

i = 10^2; %FilterLength

for iteration = 1:1

    iteration
    
    %{
    %Backward Training
    for iter1 = 1:i

            if rand-0.5 >= 0
                        x(iter1) = 1;
                    else
                        x(iter1) = -1;
            end


            yb(:,1) = H.'*( gc(:,1)*x(iter1) )+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)]; 
            vc(:,1) = vc(:,1)+StepSize*yb(:,1)*conj(x(iter1)-vc(:,1)'*yb(:,1));

    end
    %}
    
    %Normalize Transmitters
    vc(:,1)=vc(:,1)/norm(vc(:,1));



    %Forward Training
    for iter2 = 1:i
        
        iter2;

            if rand-0.5 >= 0
                        x(iter2) = 1;
                    else
                        x(iter2) = -1;
            end

            yf(:,1) = H*( vc(:,1)*x(iter2) )+ sigma*(1/sqrt(2))*[(randn(1,1)+1i*randn(1,1));(randn(1,1)+1i*randn(1,1))]; 
            gc(:,1) = gc(:,1)+StepSize*yf(:,1)*conj(x(iter2)-gc(:,1)'*yf(:,1));
            
            MSE_LMS(iter2) = 1-vc(:,1)'*H'*gc(:,1)-(vc(:,1)'*H'*gc(:,1))'+gc(:,1)'*gc(:,1)*(sigma^2)+(vc(:,1)'*H'*gc(:,1))'*(vc(:,1)'*H'*gc(:,1));
            MSE_LMS(iter2)
            SINR(iter2)= norm( gc(:,1)'*H*vc(:,1) )^2/norm( gc(:,1)'*sigma^2*gc(:,1) ); 
            C(iter2)=log2(1+SINR(iter2)); 
            MMSE(iter2) = real(   1-vc(:,1)'*H'* inv(H*vc(:,1)*vc(:,1)'*H'+eye(2)*sigma^2)*H*vc(:,1));
            
    end 
    
    SINR_w= norm( gc_w(:,1)'*H*vc(:,1) )^2/norm( gc_w(:,1)'*sigma^2*gc_w(:,1) ); 
    MMSE_w = real(   1-vc(:,1)'*H'* inv(H*vc(:,1)*vc(:,1)'*H'+eye(2)*sigma^2)*H*vc(:,1));
    %gc(:,1)=gc(:,1)/norm(gc(:,1));

        %R=H*vc(:,1)*H'*vc(:,1)'+sigma^2 
        %P=H*vc(:,1)
        %1X1
        %MSE_LMS(i) = real(   1-H'*vc(:,1)'* (1/(H*vc(:,1)*H'*vc(:,1)'+sigma^2 )) *H*vc(:,1) + ( gc(:,1) - (  1/( H*vc(:,1)*H'*vc(:,1)'+sigma^2 ) ) * H*vc(:,1)    )'*(H*vc(:,1)*H'*vc(:,1)'+1)*( gc(:,1) - (1/(H*vc(:,1)*H'*vc(:,1)'+sigma^2)) *H*vc(:,1)));
      
        

   
        %SINR_C_Wiener(i)=(   norm( G_Wiener'*H*vc(:,1) )^2   )/( norm( G_Wiener'*sigma^2*G_Wiener )    );
        %C_Wiener(i)=log2(1+SINR_C_Wiener(iter));
    
        
end
   



x=1:iter2-1;

subplot(2,1,1)
plot(x,MSE_LMS(x),x,MMSE(x))
legend('MSE(LMS)','MMSE')
xlabel('Time n')
ylabel('MSE')
title('Beam Steering;2 Antenna')

subplot(2,1,2)
plot(x,SINR(x))
legend('SINR')
xlabel('Time n')
ylabel('SINR')
title('Beam Steering;2 Antenna')


