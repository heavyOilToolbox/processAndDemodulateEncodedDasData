function GPSSignal = gpuHighBWLockInDetector(IRIGBSignal, fc, fs)
    % high bandwidth lockin detector
    B = fir1(16384, fc / fs);
    
    samplesPerCarrierPeriod = fs / fc; 
    ns = round(samplesPerCarrierPeriod / 2); % samples per half period
    dt = gpuArray((0:length(IRIGBSignal)-1) / fs); dt = dt';
    
    C0 = sin(2 * pi * fc * dt); % carrier signal
    C1 = cos(2 * pi * fc * dt); % carrier signal shifted 90 degrees out of phsae

    % high bandwidth lockin demodulator (product detector #4)
    X_shifted = gpuArray(zeros(size(IRIGBSignal)));
    X_shifted((ns + 1):end) = IRIGBSignal(1:(end - ns)); % phase shift 90 degrees forward
    Xi = IRIGBSignal .* C0;
    Xq = IRIGBSignal .* C1;
    Xis = X_shifted .* C1;
    Xqs = X_shifted .* C0;
    X1 = Xq - Xqs;
    X2 = Xi + Xis;
    
    %X1f = filtfilt(b, a, X1);
    %X2f = filtfilt(b, a, X2);
    X1f = fastFIRFilter(B, X1);
    X2f = fastFIRFilter(B, X2);
    
    GPSSignal = sqrt(2) * sqrt(X1f.^2 + X2f.^2);
    % GPSSignal = GPSSignal(padSize+1 : end-padSize);
end