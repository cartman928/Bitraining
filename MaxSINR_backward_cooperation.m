function [Vu_w, Vm_w] = MaxSINR_backward_cooperation(Z, Gu_w, Gm_w, M2, n0, Bu, Bm, upower, mpower)
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
end


%Create Big Z Matrix
for Tr_user_idx = 1:M
        for Re_user_idx = 1:M
                D{Tr_user_idx,Re_user_idx} = Z(:,:,Tr_user_idx,Re_user_idx);
        end
end   
Big_Z = cell2mat(D);

%Create Gc Vector
for user_idx = 1:M
        A{user_idx,:} = Gm_w(:,user_idx);
end  
Gc = cell2mat(A);

%Create Gp Matrix
for Tr_user_idx = 1:M
      for Re_user_idx = 1:M
          if Re_user_idx==Tr_user_idx
              B{Tr_user_idx,Re_user_idx} = Gu_w(:,Tr_user_idx)*Gu_w(:,Tr_user_idx)';
          else
              B{Tr_user_idx,Re_user_idx} = 0*eye(Nr);
          end
      end
end  
Gp = cell2mat(B);    


V = (Big_Z*Gc*Gc'*Big_Z'+Big_Z*Gp*Big_Z'+n0*eye(M*Nt))\(Big_Z*Gc);
V = sqrt(M)*V./norm(V);

for user_idx = 1 : M   
    Vm_w(:,user_idx) = V(2*user_idx-1:2*user_idx,1);
end
