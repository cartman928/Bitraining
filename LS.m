function [Gu, Gm] = LS_forward(H, Vu, Vm, M1, n0, Bu, Bm, upower, mpower)
%update filters by least square algorithm 
%notations are based on forward directions: V(transmitter),G(receiver)

[Nr,Nt,M,ignore] = size(H);

Y = zeros(Nr,M1,M);

Gu = zeros(Nr,M);
Gm = zeros(Nr,M);
for user_idx = 1 : M
    for sym_idx = 1:M1
        Y(:,sym_idx,user_idx) = H(:,:,user_idx,user_idx)*(upower(user_idx)*Vu(:,user_idx)*Bu(sym_idx,user_idx) + mpower(user_idx)*Vm(:,user_idx)*Bm(sym_idx,user_idx)) + sqrt(n0)*(randn(Nr,1)+1i*randn(Nr,1))/sqrt(2);
        
        for interf_idx = 1 : M
            if interf_idx ~= user_idx
                Y(:,sym_idx,user_idx) = Y(:,sym_idx,user_idx) +  H(:,:,user_idx,interf_idx)*(upower(interf_idx)*Vu(:,interf_idx)*Bu(sym_idx,interf_idx) + mpower(interf_idx)*Vm(:,interf_idx)*Bm(sym_idx,interf_idx));
            end
        end
    end
    
    Gu(:,user_idx) = (Y(:,:,user_idx)*Y(:,:,user_idx)')\Y(:,:,user_idx)*conj(Bu(:,user_idx));
    %Gu(:,user_idx) = Gu(:,user_idx)./norm(Gu(:,user_idx));
    
    Gm(:,user_idx) = (Y(:,:,user_idx)*Y(:,:,user_idx)')\Y(:,:,user_idx)*conj(Bm(:,user_idx));
    %Gm(:,user_idx) = Gm(:,user_idx)./norm(Gm(:,user_idx));
end

    
    