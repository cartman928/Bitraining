clc
clear
%Beamformer and Receivers for the very first iteration

vc(:,1) = [1;1]/norm([1;1]);
vp(:,1) = [1;1]/norm([1;1]);
vc(:,2) = [1;1]/norm([1;1]);
vp(:,2) = [1;1]/norm([1;1]);
gc(:,1) = [1;1]/norm([1;1]);
gp(:,1) = [1;1]/norm([1;1]);
gc(:,2) = [1;1]/norm([1;1]);
gp(:,2) = [1;1]/norm([1;1]);

%Channel
%{
H{1,1}=[0.85 0.3;0.9 0.35];
H{1,2}=[0.7 0.7;0.55 0.75];
H{2,1}=[0.9 0.5;0.73 0.6];
H{2,2}=[0.69 0.73;0.69 0.91];
%}

H{1,1}=[1 1;1 1];
H{1,2}=[1 1;1 1];
H{2,1}=[1 1;1 1];
H{2,2}=[1 1;1 1];

sigma = 10^(-2);
StepSize=10^(-5);


%{
H{1,1}=[1 1;1 1];
H{1,2}=[1 1;1 1];
H{2,1}=[1 1;1 1];
H{2,2}=[1 1;1 1];

sigma = 0.001;
%}

for iter = 1:100000 %Number of iterations

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
       
        for k = 1:2 %For user 1 and 2
            u(:,k) = [0;0];
            for j = 1:2
                    u(:,k) =  H{k,j}*( vc(:,j)*x(iter)+vp(:,j)*xp(iter,j) ) + sigma*[rand;rand];
            end
        end
        
        for k = 1:2
            
            if k==1
                m=2;
            else
                m=1;
            end
            
            gc(:,k) = gc(:,k)+StepSize*u(:,k)*(x(iter)-gc(:,k)'*u(:,k)-gc(:,m)'*u(:,m));  gc(:,k) = gc(:,k)/norm(gc(:,k))%normailize
            %gc(:,k) = gc(:,k)+StepSize*u(:,k)*(x(iter)-gc(:,k)'*u(:,k));  gc(:,k) = gc(:,k)/norm(gc(:,k));%normailize
            gp(:,k) = gp(:,k)+StepSize*u(:,k)*(xp(iter,k)-gp(:,k)'*u(:,k));  gp(:,k) = gp(:,k)/norm(gp(:,k))%normailize
            
        end
            
        %Decoding
        bc_1(iter)= sign ( gc(:,1)'*u(:,1) );
        bp_1(iter)= sign ( gp(:,1)'*u(:,1) );
        bc_2(iter)= sign ( gc(:,2)'*u(:,2) );
        bp_2(iter)= sign ( gp(:,2)'*u(:,2) );
        
        
        %Calculate Transmitters
        for k = 1:2 %For user 1 and 2
            u(:,k) = [0;0];
            for j = 1:2
                    u(:,k) =  H{j,k}.'*(gc(:,j)*x(iter)+gp(:,j)*xp(iter,j)) +sigma*[rand;rand];
            end
        end
        
        for k = 1:2
            if k == 1
                m=2;
            else
                m=1;
            end
      
            %vc(:,k) = vc(:,k)+StepSize*u(:,k)*(x(iter)-vc(:,k)'*u(:,k));
            vc(:,k) = vc(:,k)+StepSize*u(:,k)*(x(iter)-vc(:,k)'*u(:,k)-vc(:,m)'*u(:,m));
            vc(:,k) = vc(:,k)/norm(vc(:,k)); %normailize
            vp(:,k) = vp(:,k)+StepSize*u(:,k)*(xp(iter,k)-vp(:,k)'*u(:,k));
            vp(:,k) = vp(:,k)/norm(vp(:,k)); %normailize
        end
 

end

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



