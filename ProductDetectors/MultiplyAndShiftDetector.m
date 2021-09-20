function GPSSignal = MultiplyAndShiftDetector(IRIGBSignal, fc, fs, varargin)
% two-copy multiply and shift

    [IRIGBSignal, fc, fs, ~, ~, b, a] = ParseInputArguments(IRIGBSignal, fc, fs, varargin{:});

    samplesPerCarrierPeriod = fs / fc; 
    ns = round(samplesPerCarrierPeriod / 2); % samples per half period
    padSize = sqrt(fs);
    [dt, Y] = detrendAndPadSignal(IRIGBSignal, padSize, fs);
    
    C0 = sin(2 * pi * fc * dt); % carrier signal
    C1 = cos(2 * pi * fc * dt); % carrier signal shifted 90 degrees out of phsae
    S1 = C0 .* Y; % signal * carrier
    S2 = C1 .* Y; % signal * out-of-phase carrier
    S3 = zeros(size(S2)); 
    S3((ns + 1):end) = S2(1:(end - ns)); % phase-shifted copy of second product
    S1f = filtfilt(b, a, S1); 
    S3f = filtfilt(b, a, S3); 
    GPSSignal = 2*sqrt(S3f .^ 2 + S1f .^ 2); 
    GPSSignal = GPSSignal(padSize+1 : end-padSize);
end