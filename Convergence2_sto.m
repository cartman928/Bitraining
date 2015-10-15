%{

   0.1680 + 0.0371i
   0.1253 - 0.4126i

%}

clc
clear
%Beamformer and Receivers for the very first iteration

%vc(:,1) = 1;
%vc(:,1) = [1;1]/norm([1;1]);
%vp(:,1) = [1;1]/norm([1;1]);
gc(:,1) = [1;1];
gp(:,1) = [1;1];
%gc(:,1) = [1;1]/norm([1;1]);
%gp(:,1) = [1;1]/norm([1;1]);
%{
gc(:,1) = [1;1]/norm([1;1]);
gp(:,1) = [1;1]/norm([1;1]);
%}

%Channel
%{
H{1,1}=[0.85 0.3;0.9 0.35];
H{1,2}=[0.7 0.7;0.55 0.75];
H{2,1}=[0.9 0.5;0.73 0.6];
H{2,2}=[0.69 0.73;0.69 0.91];
%}

%H{1,1}=(1/sqrt(2))*[randn(1,1)+1i*randn(1,1) 0.8*(randn(1,1)+1i*randn(1,1));0.8(randn(1,1)+1i*randn(1,1)) randn(1,1)+1i*randn(1,1)];
H{1,1}=(1/sqrt(2))*[-0.9704 + 0.4012i 0.8*(1.9144 - 0.3561i);0.8*(0.4516 - 1.4800i) 0.0501 - 0.1627i];
%H{1,1}=-0.9704 + 0.4012i;

%H{1,1}=[1 0.5;0.8 1];
%H{1,1}=0.5;

sigma = 10^(-7);
StepSize=10^(-4);


%{
H{1,1}=[1 1;1 1];
H{1,2}=[1 1;1 1];
H{2,1}=[1 1;1 1];
H{2,2}=[1 1;1 1];

sigma = 0.001;
%}

    for iter = 1:10^(5) %Number of iterations

        %Random Source Message
        %Common Message
        if rand-0.5 >= 0
                    x(iter) = 1;
                else
                    x(iter) = -1;
        end

        %Private Message for User#1
        if rand-0.5 >= 0
                    xp(iter,1) = 1;
                else
                    xp(iter,1) = -1;
        end

        %Private Message for User#2
        if rand-0.5 >= 0
                    xp(iter,2) = 1;
                else
                    xp(iter,2) = -1;
        end


        %Calculate Receivers
       
        k = 1; %For user 1 and 2
            u(:,k) = [0;0];
            %u(:,k) = 0;
            j=1;
            %u(:,k) =  H{k,j}*( vc(:,j)*x(iter)+vp(:,j)*xp(iter,j) ) + sigma*[rand;rand];
            vc(:,1) = [1;1];
            vp(:,1) = [1;1];
            u(:,k) =  H{k,j}*( vc(:,j)*x(iter)+vp(:,j)*xp(iter,j) ) + sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1);randn(1,1)+1i*randn(1,1)];
            %u(:,k) =  H{k,j}*( vc(:,j)*x(iter) ) + sigma*rand;
            gc(:,k) = gc(:,k)+StepSize*u(:,k)*conj(x(iter)-gc(:,k)'*u(:,k))  
            %gc(:,k) = gc(:,k)/norm(gc(:,k));%normailize
     
            gp(:,k) = gp(:,k)+StepSize*u(:,k)*conj(xp(iter,k)-gp(:,k)'*u(:,k)) %normailize
end       
       
            
        %{
        Decoding
        bc_1(iter)= sign ( gc(:,1)'*u(:,1) );
        bp_1(iter)= sign ( gp(:,1)'*u(:,1) );
        bc_2(iter)= sign ( gc(:,2)'*u(:,2) );
        bp_2(iter)= sign ( gp(:,2)'*u(:,2) );
        %}

%{ 
for iter = 1:10^(5) %Number of iterations      
        %Calculate Transmitters
            k = 1; %For user 1 and 2
            u(:,k) = [0;0];
            %u(:,k) = 0;
            j = 1;
            %u(:,k) =  H{j,k}.'*(gc(:,j)*x(iter)+gp(:,j)*xp(iter,j)) +sigma*[rand;rand];
            
            %gc(:,j) = [1;1];
            
            u(:,k) =  H{j,k}.'*(gc(:,j)*x(iter)) +sigma*(1/sqrt(2))*[randn(1,1)+1i*randn(1,1)];
       
     
      
            %vc(:,k) = vc(:,k)+StepSize*u(:,k)*(x(iter)-vc(:,k)'*u(:,k));
            vc(:,k) = vc(:,k)+StepSize*u(:,k)*conj(x(iter)-vc(:,k)'*u(:,k))
            %vc(:,k) = vc(:,k)/norm(vc(:,k)) %normailize
            
            %vp(:,k) = vp(:,k)+StepSize*u(:,k)*(xp(iter,k)-vp(:,k)'*u(:,k));vp(:,k) = vp(:,k)/norm(vp(:,k)) %normailize
            
   
 

end
 
%}       


%{        
error=0;
for k = 101:1000
                    error = error + abs(x(k)-bc_1(k))+abs(x(k)-bc_2(k))+abs(x(k)-bc_2(k))+abs(xp(k,1)-bp_1(k))+abs(xp(k,2)-bp_2(k));
            end
4*900-error
%}

%{
x(9991:10000)
bc_1(9991:10000)
bc_2(9991:10000)

xp(9991:10000,1)'
bp_1(9991:10000)
xp(9991:10000,2)'
bp_2(9991:10000)
%}



