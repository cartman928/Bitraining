%1 user, 1X1 MIMO Channel
%calculate G,V 
%not-normalized
%calculate MSE
clc
clear

v=1.4291 - 0.6656i;
v=v/norm(v);

g=-0.9072 + 0.2915i;
g=g/norm(g);

H=-0.9704 + 0.4012i;

sigma = sqrt(10^(-3));
StepSize = 10^(-2);
xf = zeros(1,10^(7));
xb = zeros(1,10^(7));
ub = zeros(1,10^(7));
uf = zeros(1,10^(7));
SINR = zeros(1,10^(7));
C = zeros(1,10^(7));
MSE_LMS = zeros(1,10^(7));
C_LMS = zeros(1,10^(7));

for TL = 1:300
    
    TL
    
    v=1.4291 - 0.6656i;
    v=v/norm(v);

    g=-0.9072 + 0.2915i;
    g=g/norm(g);
    
   
    
    for i = 1:TL

            if rand-0.5 >= 0
                        xb(i) = 1;
                    else
                        xb(i) = -1;
            end


            ub = 0;
            ub = H*g*xb(i)+(sigma^2)*(randn(1,1)+1i*randn(1,1));
             v = v+StepSize*ub*conj(xb(i)-v'*ub);

    end
    
    for i = 1:TL

            if rand-0.5 >= 0
                        xf(i) = 1;
                    else
                        xf(i) = -1;
            end

            uf = 0;
            uf = H*v*xf(i)+(sigma^2)*(randn(1,1)+1i*randn(1,1));
             g = g+StepSize*uf*conj(xf(i)-g'*uf);

    end
    
    
    SINR(TL)=  norm( g'*H*v )^2/norm( g'*sigma^2*g );
    C(TL)=log2(1+SINR(TL));
      
end

    x=1:TL;
 
    plot(x,C(x))
    legend('C')
    xlabel('Iteration')
    ylabel('C(bit/channel)')
    title('1User;1X1 SISO;Only Common Msg;M=1')



