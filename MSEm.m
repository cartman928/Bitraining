function E = MSEm(H, Gu, Gm, Vu, Vm, n0, upower, mpower)
%update receive filters by least square algorithm

[Nr,Nt,M,ignore] = size(H);

Y = zeros(1,Nr);


for user_idx = 1 : M
        Y(user_idx) = 1-mpower(user_idx)*Vm(:,user_idx)'*H(:,:,user_idx,user_idx)'*mpower(user_idx)*Gm(:,user_idx)-(mpower(user_idx)*Vm(:,user_idx)'*H(:,:,user_idx,user_idx)'*mpower(user_idx)*Gm(:,user_idx))'...
                       +mpower(user_idx)*Gm(:,user_idx)'*H(:,:,user_idx,user_idx)*mpower(user_idx)*Vm(:,user_idx)*((user_idx)*Gm(:,user_idx)'*H(:,:,user_idx,user_idx)*mpower(user_idx)*Vm(:,user_idx))'...
                       +mpower(user_idx)*Gm(:,user_idx)'*H(:,:,user_idx,user_idx)*upower(user_idx)*Vu(:,user_idx)*(mpower(user_idx)*Gm(:,user_idx)'*H(:,:,user_idx,user_idx)*upower(user_idx)*Vu(:,user_idx))'...
                       +mpower(user_idx)*Gm(:,user_idx)'*n0*eye(Nr)*mpower(user_idx)*Gm(:,user_idx);
        
        interferencem = zeros(Nt,1);
        for interf_idx = 1 : M
            if interf_idx ~= user_idx
                
                interferencem = interferencem + H(:,:,user_idx,interf_idx)*mpower(interf_idx)*Vm(:,interf_idx);
                Y(user_idx) = Y(user_idx) -  mpower(interf_idx)*Vm(:,interf_idx)'*H(:,:,user_idx,interf_idx)'*mpower(user_idx)*Gm(:,user_idx) - mpower(user_idx)*Gm(:,user_idx)'*H(:,:,user_idx,interf_idx)*mpower(interf_idx)*Vm(:,interf_idx)...
                                          +  mpower(user_idx)*Gm(:,user_idx)'*H(:,:,user_idx,interf_idx)*upower(interf_idx)*Vu(:,interf_idx)*upower(interf_idx)*Vu(:,interf_idx)'*H(:,:,user_idx,interf_idx)'*mpower(user_idx)*Gm(:,user_idx)...
                                          +  mpower(user_idx)*Gm(:,user_idx)'*H(:,:,user_idx,user_idx)*mpower(user_idx)*Vm(:,user_idx)*mpower(interf_idx)*Vm(:,interf_idx)'*H(:,:,user_idx,interf_idx)'*mpower(user_idx)*Gm(:,user_idx)...
                                          +  mpower(user_idx)*Gm(:,user_idx)'*H(:,:,user_idx,interf_idx)*mpower(interf_idx)*Vm(:,interf_idx)*mpower(user_idx)*Vm(:,user_idx)'*H(:,:,user_idx,user_idx)'*mpower(user_idx)*Gm(:,user_idx);
                                                  
                
            end
        end

    Y(user_idx) = abs(  Y(user_idx) + mpower(user_idx)*Gm(:,user_idx)'*(interferencem*interferencem')*mpower(user_idx)*Gm(:,user_idx) );
    
end

E = sum(Y);

