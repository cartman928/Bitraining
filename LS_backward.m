function [Vu, Vm] = LS_backward(H, Gu, Gm, M2, n0, Bu, Bm, upower, mpower)
%update receive filters by least square algorithm

[Nr,Nt,M,ignore] = size(H);

Y = zeros(Nt,M2,M);

Vu = zeros(Nt,M);
Vm = zeros(Nt,M);
for user_idx = 1 : M
    for sym_idx = 1:M2
        Y(:,sym_idx,user_idx) = H(:,:,user_idx,user_idx).'*(upower(user_idx)*conj(Gu(:,user_idx))*Bu(sym_idx,user_idx) + mpower(user_idx)*conj(Gm(:,user_idx))*Bm(sym_idx,user_idx)) + sqrt(n0)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
        for interf_idx = 1 : M
            if interf_idx ~= user_idx
                Y(:,sym_idx,user_idx) = Y(:,sym_idx,user_idx) +  H(:,:,interf_idx, user_idx).'*(upower(interf_idx)*conj(Gu(:,interf_idx))*Bu(sym_idx,interf_idx) + mpower(interf_idx)*conj(Gm(:,interf_idx))*Bm(sym_idx,interf_idx));
            end
        end
    end
    
    Vu(:,user_idx) = (Y(:,:,user_idx)*Y(:,:,user_idx)')\Y(:,:,user_idx)*conj(Bu(:,user_idx));
    Vu(:,user_idx) = conj(Vu(:,user_idx));
    Vu(:,user_idx) = Vu(:,user_idx)./norm(Vu(:,user_idx));
    
    Vm(:,user_idx) = inv(Y(:,:,user_idx)*Y(:,:,user_idx)')*Y(:,:,user_idx)*conj(Bm(:,user_idx));
    Vm(:,user_idx) = conj(Vm(:,user_idx));
    Vm(:,user_idx) = Vm(:,user_idx)./norm(Vm(:,user_idx));
    
    
%     if (numiters == 4 && user_idx == 1)
%                     mean(abs(Bu(:,1).' - Vu(:,1).'*Y(:,:,1)))
%                 end
end