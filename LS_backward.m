function [Vu, Vm] = LS_backward(H, Gu, Gm, M2, n0, Bu, Bm, upower, mpower)
%update receive filters by least square algorithm
%[Vu, Vm] = LS_backward(H, Gu, Gm, M2, n0, Bbw, BbwBr, upower, mpower);


[Nr,Nt,M,ignore] = size(H);

%M=2, M2=10

for user_idx = 1 : M
    for sym_idx = 1:M2
        Y(:,sym_idx,user_idx) = H(:,:,user_idx,user_idx)'*Gu(:,user_idx)*Bbw(sym_idx,user_idx) + sqrt(10^(-3))*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
        for interf_idx = 1 : M
            if interf_idx ~= user_idx
                Y(:,sym_idx,user_idx) = Y(:,sym_idx,user_idx) +  H(:,:,interf_idx, user_idx)'*Gu(:,interf_idx)*Bbw(sym_idx,interf_idx);
            end
        end
    end    
    Vu(:,user_idx) = inv(Y(:,:,user_idx)*Y(:,:,user_idx)')*Y(:,:,user_idx)*conj(Bbw(:,user_idx));
    Vu(:,user_idx) = Vu(:,user_idx)./norm(Vu(:,user_idx));
    Vm(:,user_idx) = [0;0];
end