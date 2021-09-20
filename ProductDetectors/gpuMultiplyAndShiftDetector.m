function GPSSignal = gpuMultiplyAndShiftDetector(IRIGBSignal, fc, fs)
% two-copy multiply and shift

    samplesPerCarrierPeriod = fs / fc; 
    ns = round(samplesPerCarrierPeriod / 2); % samples per half period
    dt = gpuArray((0:length(IRIGBSignal)-1) / fs); dt = dt';
    
    C0 = sin(2 * pi * fc * dt); % carrier signal
    C1 = cos(2 * pi * fc * dt); % carrier signal shifted 90 degrees out of phsae
    S1 = C0 .* IRIGBSignal; % signal * carrier
    S2 = C1 .* IRIGBSignal; % signal * out-of-phase carrier
    S3 = gpuArray(zeros(size(S2))); 
    S3((ns + 1):end) = S2(1:(end - ns)); % phase-shifted copy of second product
    
    B = fir1(16384, fc / fs);
    %S1f = filtfilt(b, a, S1); 
    %S3f = filtfilt(b, a, S3); 
    S1f = fastFIRFilter(B, S1);
    S3f = fastFIRFilter(B, S3);
    
    GPSSignal = 2*sqrt(S3f .^ 2 + S1f .^ 2); 
    % GPSSignal = GPSSignal(padSize+1 : end-padSize);
end