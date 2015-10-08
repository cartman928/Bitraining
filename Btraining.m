%{
Communication Systems: 
Source Messages -> Beamformer -> Channel -> Receiver-> Decoded Messages
%}

%Source Message

x = 1; %Common Message
xp_1=1; %Private Message for User#1
xp_2=1; %Private Message for User#2

%Beamformer

vc_1 = [1;1];
vp_1 = [1;1];
vc_2 = [1;1];
vp_2 = [1;1];

%Signals After Beamformers

TX1 = x.*vc_1+xp_1.*vp_1;
TX2 = x.*vc_2+xp_2.*vp_2;

%Channel

H11=1;
H12=1;
H21=1;
H22=1;
H = [H11 H12; H21 H22];

%Signals Before Receivers 

dummy_RX1 = (H*TX1)+(H*TX2);
RX1 = dummy_RX1(1,1);
dummy_RX2 = (H*TX1)+(H*TX2);
RX2 = dummy_RX2(2,1);

%Receivers

gc_1 = [1;1];
gp_1 = [1;1];
gc_2 = [1;1];
gp_2 = [1;1];

%Decoded Messages

x1 = gc_1*RX1; 
xp_1=gp_1*RX1; 
x2 = gc_2*RX2; 
xp_2=gp_2*RX2; 





