function [Gu, Gm] = LS_forward(H, Vu, Vm, M1, n0, Bu, Bm, upower, mpower)
%update receive filters by least square algorithm

[Nr,Nt,M,ignore] = size(H);


for user_idx = 1 : M
    for sym_idx = 1:M1
        Y(:,sym_idx,user_idx) = H(:,:,user_idx,user_idx)*(Vu(:,user_idx)*Bbw(sym_idx,user_idx)  ) + sqrt(n0)*(randn(Nr,1)+1i*randn(Nr,1))/sqrt(2);
        for interf_idx = 1 : M
            if interf_idx ~= user_idx
                Y(:,sym_idx,user_idx) = Y(:,sym_idx,user_idx) +  H(:,:,user_idx,interf_idx)*Vu(:,interf_idx)*Bbw(sym_idx,interf_idx) ;
            end
        end
        
    end
    
    Gu(:,user_idx) = inv(Y(:,:,user_idx)*Y(:,:,user_idx)')*Y(:,:,user_idx)*conj(Bbw(:,user_idx));
    Gu(:,user_idx) = Gu(:,user_idx)./norm(Gu(:,user_idx));
    Gm(:,user_idx) = [0;0];
end

    
    