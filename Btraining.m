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

gc_1 = [1;1];
gp_1 = [1;1];
gc_2 = [1;1];
gp_2 = [1;1];

sigma = 0.2;
k=1;


sum0(:,k) = [0;0];
for j = 1:2
    if j~=k
        sum0(:,k) = sum0(:,k) + H(k,j)*vc(:,j);
    end
end

gc(:,k)=inv(  H(k,k)*vc(:,k)*vc(:,k)'*H(k,k)'+sum0(:,k)*sum0(:,k)'+H(k,k)*vc(:,k)*sum0(:,k)'+sum0(:,k)*)*vc(:,k)'*H(k,k)'+eye(2)*sigma^2)





%k=2;
%gc(k)=2;


%Decoded Messages

x1 = gc_1*RX1; 
xp_1=gp_1*RX1; 
x2 = gc_2*RX2; 
xp_2=gp_2*RX2;





