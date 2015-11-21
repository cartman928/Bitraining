function [Gu_w, Gm_w] = MaxSINR(H, Vu_w, Vm_w, n0, upower, mpower)
%update filters by MaxSINR algorithm
%notations are based on forward directions: V(transmitter),G(receiver)

[Nr,Nt,M,ignore] = size(H);

Y = zeros(Nr,Nr);

Gu_w = zeros(Nr,M);
Gm_w = zeros(Nr,M);
for user_idx = 1 : M
        signal = H(:,:,user_idx,user_idx)*mpower(user_idx)*Vm_w(:,user_idx);
        Y(:,:,user_idx) = H(:,:,user_idx,user_idx)*(mpower(user_idx)^2*Vm_w(:,user_idx)*Vm_w(:,user_idx)' + upower(user_idx)^2*Vu_w(:,user_idx)*Vu_w(:,user_idx)')*H(:,:,user_idx,user_idx)' + n0*eye(Nr);
        
        interferencem = zeros(Nr,1);
        for interf_idx = 1 : M
            if interf_idx ~= user_idx
                
                signal = signal + H(:,:,user_idx,interf_idx)*mpower(interf_idx)*Vm_w(:,interf_idx);
                interferencem = interferencem + H(:,:,user_idx,interf_idx)*mpower(interf_idx)*Vm_w(:,interf_idx);
                Y(:,:,user_idx) = Y(:,:,user_idx) +  H(:,:,user_idx,interf_idx)*(upower(interf_idx)^2*Vu_w(:,interf_idx)*Vu_w(:,interf_idx)')*H(:,:,user_idx,interf_idx)'...
                                                  +  H(:,:,user_idx,user_idx)*mpower(user_idx)*Vm_w(:,user_idx)*mpower(interf_idx)*Vm_w(:,interf_idx)'*H(:,:,user_idx,interf_idx)'...
                                                  +  (H(:,:,user_idx,user_idx)*mpower(user_idx)*Vm_w(:,user_idx)*mpower(interf_idx)*Vm_w(:,interf_idx)'*H(:,:,user_idx,interf_idx)')';
                
            end
        end

    Y(:,:,user_idx) = Y(:,:,user_idx) + interferencem*interferencem';
    
    Gu_w(:,user_idx) = Y(:,:,user_idx)\(H(:,:,user_idx,user_idx)*upower(user_idx)*Vu_w(:,user_idx));
    if norm(Gu_w(:,user_idx)) ~= 0; %prevent divided by zero
    %Gu_w(:,user_idx) = Gu_w(:,user_idx)./norm(Gu_w(:,user_idx));
    end
    
    Gm_w(:,user_idx) = Y(:,:,user_idx)\signal;
    if norm(Gm_w(:,user_idx)) ~= 0; %prevent divided by zero
    %Gm_w(:,user_idx) = Gm_w(:,user_idx)./norm(Gm_w(:,user_idx));
    end
end

    
    