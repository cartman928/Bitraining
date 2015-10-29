%1 user, 4X4 MIMO Channel
%calculate V with cooperation
%Only Common Message
%It won't converge to wiener solution. However, if we turn on private channel, it converges.

clc
clear


V = [1.4291 - 0.6656i;0.0183 - 1.1054i;0.8221 - 0.3290i;0.3581 + 0.1167i];
V=sqrt(2)*(V/norm(V));

gc(:,1) = [-0.9072 + 0.2915i;-1.1637 - 0.6096i];
gc(:,2) = [1.0489 - 0.1224i;-1.7782 - 2.3879i];
gc(:,1) = gc(:,1)/norm(gc(:,1));
gc(:,2) = gc(:,2)/norm(gc(:,2));

Gc = [gc(:,1);gc(:,2)];



K{1,1}=(1/sqrt(2))*[-0.9704 + 0.4012i (1.9144 - 0.3561i);(0.4516 - 1.4800i) 0.0501 - 0.1627i];
K{1,2}=0.8*(1/sqrt(2))*[-0.6337 + 0.8001i (-1.1485 + 0.2689i);(0.4516 - 1.4800i) -0.6651 - 0.9268i];
K{2,1}=0.8*(1/sqrt(2))*[-0.0165 + 0.9251i (-1.0463 - 0.5763i);(-0.1497 - 1.5829i) -0.7804 - 0.4109i];
K{2,2}=(1/sqrt(2))*[-0.9313 + 0.8060i (0.7313 + 0.0698i);(-0.2850 + 1.1345i) -0.5113 - 0.1662i];

H = [K{1,1} K{1,2};K{2,1} K{2,2}];


sigma = sqrt(10^(-3));
StepSize = 10^(-7);

U = [0;0;0;0];
x = zeros(1,10^(7));

for iter = 1:10^(7) 
        if rand-0.5 >= 0
                    x(iter) = 1;
                else
                    x(iter) = -1;
        end

            
        U =  H*Gc*x(iter)+ sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
         
        V = V+StepSize*U*conj(x(iter)-(V)'*U)
        
        

end

V = sqrt(2)*V/norm(V);


V_wiener = inv( H*Gc*(Gc)'*(H)'+eye(4)*sigma^2 )*H*Gc;
V_wiener = sqrt(2)*V_wiener/norm(V_wiener);



SINR_C(1)=(   norm( gc(:,1)'*K{1,1}*[V(1,1);V(2,1)]+gc(:,2)'*K{1,2}*[V(3,1);V(4,1)] )^2   )/( norm( gc(:,1)'*sigma^2*gc(:,1) )    );
SINR_C(2)=(   norm( gc(:,2)'*K{2,1}*[V(1,1);V(2,1)]+gc(:,2)'*K{2,2}*[V(3,1);V(4,1)] )^2   )/( norm( gc(:,2)'*sigma^2*gc(:,2) )    );
C_LMS=log2(1+SINR_C(1))+log2(1+SINR_C(2)) 

SINR_C_Wiener(1)=(   norm( gc(:,1)'*K{1,1}*[V_wiener(1,1);V_wiener(2,1)]+gc(:,2)'*K{1,2}*[V_wiener(3,1);V_wiener(4,1)] )^2   )/( norm( gc(:,1)'*sigma^2*gc(:,1) )    );
SINR_C_Wiener(2)=(   norm( gc(:,2)'*K{2,1}*[V_wiener(1,1);V_wiener(2,1)]+gc(:,2)'*K{2,2}*[V_wiener(3,1);V_wiener(4,1)] )^2   )/( norm( gc(:,2)'*sigma^2*gc(:,2) )    );
C_Wiener=log2(1+SINR_C_Wiener(1))+log2(1+SINR_C_Wiener(2)) 




%{

 V_wiener =

  -0.0997 - 0.0113i
  -0.0540 + 0.2831i
  -0.0669 - 0.0119i
   0.1457 + 0.4039i

%}


