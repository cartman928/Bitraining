function [Gu_w, Gm_w] = Wiener_forward(H, Vu_w, Vm_w, M1, n0, Bu, Bm, upower, mpower)
%update receive filters by least square algorithm

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
                interferencem = interferencem + H(:,:,user_idx,interf_idx)*upower(interf_idx)*Vu_w(:,interf_idx);
                Y(:,:,user_idx) = Y(:,:,user_idx) +  H(:,:,user_idx,interf_idx)*(upower(interf_idx)^2*Vu_w(:,interf_idx)*Vu_w(:,interf_idx)')*H(:,:,user_idx,interf_idx)'...
                                                  +  H(:,:,user_idx,user_idx)*mpower(user_idx)*Vm_w(:,user_idx)*mpower(interf_idx)*Vm_w(:,interf_idx)'*H(:,:,user_idx,interf_idx)'...
                                                  +  (H(:,:,user_idx,user_idx)*mpower(user_idx)*Vm_w(:,user_idx)*mpower(interf_idx)*Vm_w(:,interf_idx)'*H(:,:,user_idx,interf_idx)')';
                
            end
        end

    Y(:,:,user_idx) = Y(:,:,user_idx) + interferencem*interferencem';
    
    Gu_w(:,user_idx) = Y(:,:,user_idx)\(H(:,:,user_idx,user_idx)*upower(user_idx)*Vu_w(:,user_idx));
    Gu_w(:,user_idx) = Gu_w(:,user_idx)./norm(Gu_w(:,user_idx));
    
    Gm_w(:,user_idx) = Y(:,:,user_idx)\signal;
    Gm_w(:,user_idx) = Gm_w(:,user_idx)./norm(Gm_w(:,user_idx));
end

    
    