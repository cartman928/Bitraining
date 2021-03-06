function [Vu, Vm] = LS_cooperation(Z, Gu, Gm, M2, n0, Bu, Bm, upower, mpower)
%update filters by least square algorithm
%notations are based on backward directions: G(transmitter),V(receiver)
[Nr,Nt,M,ignore] = size(Z);

Y = zeros(Nt,M2,M);

Vu = zeros(Nt,M);
Vm = zeros(Nt,M);
for user_idx = 1 : M
    for sym_idx = 1:M2
        Y(:,sym_idx,user_idx) = Z(:,:,user_idx,user_idx)*(upower(user_idx)*Gu(:,user_idx)*Bu(sym_idx,user_idx) + mpower(user_idx)*Gm(:,user_idx)*Bm(sym_idx,user_idx)) + sqrt(n0)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
        for interf_idx = 1 : M
            if interf_idx ~= user_idx
                Y(:,sym_idx,user_idx) = Y(:,sym_idx,user_idx) +  Z(:,:,user_idx,interf_idx)*(upower(interf_idx)*Gu(:,interf_idx)*Bu(sym_idx,interf_idx) + mpower(interf_idx)*Gm(:,interf_idx)*Bm(sym_idx,interf_idx));
            end
        end
    end
    
    Vu(:,user_idx) = (Y(:,:,user_idx)*Y(:,:,user_idx)')\Y(:,:,user_idx)*conj(Bu(:,user_idx));
    Vu(:,user_idx) = Vu(:,user_idx)./norm(Vu(:,user_idx));
    
    %Create Shared Data Matrix
    D{user_idx,1} = Y(:,:,user_idx); 
       
end


Big_Y = cell2mat(D);
V = (Big_Y*Big_Y')\Big_Y*conj(Bm(:,1));
if norm(V) ~= 0;
%V = sqrt(M)*V./norm(V);
end

for user_idx = 1 : M   
    Vm(:,user_idx) = V(2*user_idx-1:2*user_idx);
end