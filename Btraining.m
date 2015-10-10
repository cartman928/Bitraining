%{
Communication Systems: 
Source Messages -> Beamformer -> Channel -> Receiver-> Decoded Messages
%}

%Source Message

x = 1; %Common Message
xp_1=1; %Private Message for User#1
xp_2=1; %Private Message for User#2

%Beamformer

vc(:,1) = [1;1];
vp(:,1) = [1;1];
vc(:,2) = [1;1];
vp(:,2) = [1;1];

%Signals After Beamformers

TX(:,1) = x*vc(:,1)+xp_1*vp(:,1);
TX(:,2) = x*vc(:,2)+xp_2*vp(:,2);

%Channel

H(1,1)=1;
H(1,2)=1;
H(2,1)=1;
H(2,2)=1;

%Signals Before Receivers 
noise = 1;
dummy_RX1 = (H*TX(:,1))+(H*TX(:,2));
RX1 = dummy_RX1(1,1)+noise;
dummy_RX2 = (H*TX(:,1))+(H*TX(:,2));
RX2 = dummy_RX2(2,1)+noise;

%Receivers


sigma = 0.2;%Square Root of Noise Variance
k=1;


sum_c1(:,k) = [0;0];
for j = 1:2
    if j~=k
        sum_c1(:,k) = sum_c1(:,k) + H(k,j)*vc(:,j);
    end
end

sum_p1(:,k) = [0;0];
for j = 1:2
    if j~=k
        sum_p1(:,k) = sum_p1(:,k) + H(k,j)*vp(:,j);
    end
end

sum_c2(:,k) = [0;0];
for j = 1:2
        sum_c2(:,k) = sum_c2(:,k) + H(k,j)*vc(:,j);
end




gc(:,k) = inv(  H(k,k)*vc(:,k)*vc(:,k)'*H(k,k)'+sum_c1(:,k)*sum_c1(:,k)'+H(k,k)*vc(:,k)*sum_c1(:,k)'+sum_c1(:,k)*vc(:,k)'*H(k,k)'+eye(2)*sigma^2+H(k,k)*vp(:,k)*vp(:,k)'*H(k,k)'+sum_p1(:,k)*sum_p1(:,k)'  ) * ( sum_c2(:,k) );
gp(:,k) = inv(  H(k,k)*vc(:,k)*vc(:,k)'*H(k,k)'+sum_c1(:,k)*sum_c1(:,k)'+H(k,k)*vc(:,k)*sum_c1(:,k)'+sum_c1(:,k)*vc(:,k)'*H(k,k)'+eye(2)*sigma^2+H(k,k)*vp(:,k)*vp(:,k)'*H(k,k)'+sum_p1(:,k)*sum_p1(:,k)'  ) * ( H(k,k)*vp(:,k) );

vc(:,k) = -2*gc(:,k)-2*lambda(k)



vp(:,k)=;


 %k=2;
%gc(k)=2;


%Decoded Messages

x1 = gc_1*RX1; 
xp_1=gp_1*RX1; 
x2 = gc_2*RX2; 
xp_2=gp_2*RX2;





