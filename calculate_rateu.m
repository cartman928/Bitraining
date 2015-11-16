function sumrate = calculate_rateu (h, n0, beamformeru, filteru, beamformerm, upower, mpower)
%calculate the sum rate

[Nr,Nt,Nu,ignore] = size(h);
rate = zeros(1, Nu);


for k = 1:Nu
    beamformeru(:,k) = beamformeru(:,k)./norm(beamformeru(:,k));
    filteru(:, k) = filteru(:, k)./norm(filteru(:,k));
end


for k = 1:Nu
    signal = upower(k).^2*abs(filteru(:,k)'* h(:,:,k,k)*beamformeru(:,k))^2;
    if signal < eps
        rate(k) = 0;
        continue
    end
    interference = 0;
    interferencem = mpower(k)*filteru(:,k)'* h(:,:,k,k)*beamformerm(:,k);
    for interf_idx = 1 : Nu
        if interf_idx ~= k
            interference = interference + abs(filteru(:,k)'* h(:,:,k,interf_idx)*(upower(interf_idx)*beamformeru(:,interf_idx)))^2;
            interferencem = interferencem + filteru(:,k)'* h(:,:,k,interf_idx)*(mpower(interf_idx)*beamformerm(:,interf_idx));
        end
    end
    interference = interference + abs(interferencem)^2;
    noise = filteru(:,k)'*filteru(:,k)*n0;
    rate(k) = log2(1+signal/(noise+interference));
    
    
end



rate;

sumrate = sum(rate);

    