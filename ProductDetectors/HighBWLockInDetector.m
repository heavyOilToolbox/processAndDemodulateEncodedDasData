function GPSSignal = HighBWLockInDetector(IRIGBSignal, fc, fs, varargin)
    % high bandwidth lockin detector
    [IRIGBSignal, fc, fs, ~, ~, b, a] = ParseInputArguments(IRIGBSignal, fc, fs, varargin{:});
    samplesPerCarrierPeriod = fs / fc; 
    ns = round(samplesPerCarrierPeriod / 2); % samples per half period
    padSize = sqrt(fs);
    [dt, Y] = detrendAndPadSignal(IRIGBSignal, padSize, fs);
    C0 = sin(2 * pi * fc * dt); % carrier signal
    C1 = cos(2 * pi * fc * dt); % carrier signal shifted 90 degrees out of phsae

    % high bandwidth lockin demodulator (product detector #4)
    X_shifted = zeros(size(Y));
    X_shifted((ns + 1):end) = Y(1:(end - ns)); % phase shift 90 degrees forward
    Xi = Y .* C0;
    Xq = Y .* C1;
    Xis = X_shifted .* C1;
    Xqs = X_shifted .* C0;
    X1 = Xq - Xqs;
    X2 = Xi + Xis;
    X1f = filtfilt(b, a, X1);
    X2f = filtfilt(b, a, X2);
    GPSSignal = sqrt(2) * sqrt(X1f.^2 + X2f.^2);
    GPSSignal = GPSSignal(padSize+1 : end-padSize);
end