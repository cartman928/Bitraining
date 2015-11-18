function [Vu_w, Vm_w] = MaxSINR_backward(Z, Gu_w, Gm_w, n0, upower, mpower)
%update receive filters by least square algorithm

[Nr,Nt,M,ignore] = size(Z);

Y = zeros(Nt,Nt);

Vu_w = zeros(Nr,M);
Vm_w = zeros(Nr,M);
for user_idx = 1 : M
        signal = Z(:,:,user_idx,user_idx)*mpower(user_idx)*Gm_w(:,user_idx);
        Y(:,:,user_idx) = Z(:,:,user_idx,user_idx)*(    mpower(user_idx)^2*Gm_w(:,user_idx)*Gm_w(:,user_idx)' + upower(user_idx)^2*Gu_w(:,user_idx)*Gu_w(:,user_idx)'   )*Z(:,:,user_idx,user_idx)' + n0*eye(Nr);
        
        interferencem = zeros(Nt,1);
        for interf_idx = 1 : M
            if interf_idx ~= user_idx
                
                signal = signal + Z(:,:,user_idx,interf_idx)*mpower(interf_idx)*Gm_w(:,interf_idx);
                interferencem = interferencem + Z(:,:,user_idx,interf_idx)*mpower(interf_idx)*Gm_w(:,interf_idx);
                Y(:,:,user_idx) = Y(:,:,user_idx) +  Z(:,:,user_idx,interf_idx)*(upower(interf_idx)^2*Gu_w(:,interf_idx)*Gu_w(:,interf_idx)')*Z(:,:,user_idx,interf_idx)'...
                                                  +  Z(:,:,user_idx,user_idx)*mpower(user_idx)*Gm_w(:,user_idx)*mpower(interf_idx)*Gm_w(:,interf_idx)'*Z(:,:,user_idx,interf_idx)'...
                                                  +  (Z(:,:,user_idx,user_idx)*mpower(user_idx)*Gm_w(:,user_idx)*mpower(interf_idx)*Gm_w(:,interf_idx)'*Z(:,:,user_idx,interf_idx)')';
                
            end
        end

    Y(:,:,user_idx) = Y(:,:,user_idx) + interferencem*interferencem';
    
    Vu_w(:,user_idx) = Y(:,:,user_idx)\(Z(:,:,user_idx,user_idx)*upower(user_idx)*Gu_w(:,user_idx));
    Vu_w(:,user_idx) = Vu_w(:,user_idx)./norm(Vu_w(:,user_idx));
    
    Vm_w(:,user_idx) = Y(:,:,user_idx)\signal;
    Vm_w(:,user_idx) = Vm_w(:,user_idx)./norm(Vm_w(:,user_idx));
end

    
    