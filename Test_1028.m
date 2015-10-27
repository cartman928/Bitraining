%2 user, 2X2 MIMO Channel
%calculate V with cooperation
%Turn off priviate channel
%It won't converge to wiener solution. However, if we turn on private channel, it converges.

clc
clear


%V = [1.4291 - 0.6656i;0.0183 - 1.1054i;0.8221 - 0.3290i;0.3581 + 0.1167i];
V = [1.4291 - 0.6656i;0.0183 - 1.1054i];
%V=sqrt(2)*(V/norm(V));
V=(V/norm(V));

%Gc = [-0.9072 + 0.2915i;-1.1637 - 0.6096i;1.0489 - 0.1224i;-1.7782 - 2.3879i];
Gc = [-0.9072 + 0.2915i;-1.1637 - 0.6096i];
%Gc=sqrt(2)*(Gc/norm(Gc));
Gc=(Gc/norm(Gc));

K{1,1}=(1/sqrt(2))*[-0.9704 + 0.4012i (1.9144 - 0.3561i);(0.4516 - 1.4800i) 0.0501 - 0.1627i];
K{1,2}=0.8*(1/sqrt(2))*[-0.6337 + 0.8001i (-1.1485 + 0.2689i);(0.4516 - 1.4800i) -0.6651 - 0.9268i];
K{2,1}=0.8*(1/sqrt(2))*[-0.0165 + 0.9251i (-1.0463 - 0.5763i);(-0.1497 - 1.5829i) -0.7804 - 0.4109i];
K{2,2}=(1/sqrt(2))*[-0.9313 + 0.8060i (0.7313 + 0.0698i);(-0.2850 + 1.1345i) -0.5113 - 0.1662i];

%H = [K{1,1} K{1,2};K{2,1} K{2,2}];
H = [K{1,1}];

sigma = sqrt(10^(-3));
StepSize = 10^(-6);

U = [0;0;0;0];
x = zeros(1,2*10^(6));

for iter = 1:2*10^(6) 
        if rand-0.5 >= 0
                    x(iter) = 1;
                else
                    x(iter) = -1;
        end

            
            U =  H*Gc*x(iter)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
         
       V = V+StepSize*U*conj(x(iter)-(V)'*U)

end

V_wiener = inv( H*Gc*(Gc)'*(H)'+eye(2)*sigma^2 )*H*Gc;

%{

 V_wiener =

  -0.5162 - 0.4296i
  -0.0417 + 0.5013i
%}


