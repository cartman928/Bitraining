function [Gu_w, Gm_w] = MaxSINRm_only(H, Vu_w, Vm_w, n0, upower, mpower)
%update filters by MaxSINR algorithm
%notations are based on forward directions: V(transmitter),G(receiver)

[Nr,Nt,M,ignore] = size(H);

Y = zeros(Nr,Nr);

Gu_w = zeros(Nr,M);
Gm_w = zeros(Nr,M);
for user_idx = 1 : M
        for i = 1 : M
        Y(:,user_idx) = Y(:,user_idx)+H(:,:,user_idx,i)*Vm_w(:,i);
        end
    
    Gu_w(:,user_idx) = zeros(Nr,1);
    if norm(Gu_w(:,user_idx)) ~= 0; %prevent divided by zero
    %Gu_w(:,user_idx) = Gu_w(:,user_idx)./norm(Gu_w(:,user_idx));
    end
    
    Gm_w(:,user_idx) = (Y(:,user_idx)*Y(:,user_idx)'+eye(Nr)*n0)\Y(:,user_idx);
    if norm(Gm_w(:,user_idx)) ~= 0; %prevent divided by zero
    %Gm_w(:,user_idx) = Gm_w(:,user_idx)./norm(Gm_w(:,user_idx));
    end
end

    
    