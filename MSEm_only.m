function E = MSEm_only(H, Gu, Gm, Vu, Vm, n0, upower, mpower)
%update receive filters by least square algorithm
%notations are based on forward directions: V(transmitter),G(receiver)

[Nr,Nt,M,ignore] = size(H);

Y = zeros(1,Nr);


for user_idx = 1 : M
        Y(user_idx) = 1 +Gm(:,user_idx)'*Gm(:,user_idx)*n0;
        
        signal = zeros(Nt,1);
        for i = 1 : M       
                signal = signal + H(:,:,user_idx,i)*Vm(:,i);
                Y(user_idx) = Y(user_idx) -  Vm(:,i)'*H(:,:,user_idx,i)'*Gm(:,user_idx) - (Vm(:,i)'*H(:,:,user_idx,i)'*Gm(:,user_idx))';
        end

    Y(user_idx) = abs(  Y(user_idx) + Gm(:,user_idx)'*(signal*signal')*Gm(:,user_idx) );
    
end

E = sum(Y);

