function sumrate = calculate_ratem (h, n0, beamformerm, filterm, beamformeru, upower, mpower)
%calculate the sum rate

[Nr,Nt,Nu,ignore] = size(h);
rate = zeros(1, Nu);


for k = 1:Nu
    beamformerm(:,k) = beamformerm(:,k)./norm(beamformerm(:,k));
    filterm(:, k) = filterm(:, k)./norm(filterm(:,k));
end


for k = 1:Nu
    signal = mpower(k)*filterm(:,k)'* h(:,:,k,k)*beamformerm(:,k);
    
    
    
    
    interference = 0;
    interference = upower(k).^2*abs(filterm(:,k)'* h(:,:,k,k)*beamformeru(:,k))^2;
    for interf_idx = 1 : Nu
        if interf_idx ~= k
            signal = signal + filterm(:,k)'* h(:,:,k,interf_idx)*(mpower(interf_idx)*beamformerm(:,interf_idx));
            interference = interference + abs(filterm(:,k)'* h(:,:,k,interf_idx)*(upower(interf_idx)*beamformeru(:,interf_idx)))^2;
        end
    end
    signal = abs(signal)^2;
    if signal < eps
        rate(k) = 0;
        continue
    end
    noise = filterm(:,k)'*filterm(:,k)*n0;
    rate(k) = log2(1+signal/(noise+interference));
    
    
end

rate;

sumrate = sum(rate);

