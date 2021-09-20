function GPSSignal = BasicProductDetector(IRIGBSignal, fc, fs, varargin)
    % basic product detector (product detector #3)
    [IRIGBSignal, fc, fs, INI_PHASE, ~, b, a] = ParseInputArguments(IRIGBSignal, fc, fs, varargin{:});
    initialPhase = pi * INI_PHASE;
    padSize = sqrt(fs);
    [dt, Y] = detrendAndPadSignal(IRIGBSignal, padSize, fs);
    C2 = sin(2 * pi * fc * dt + initialPhase); % carrier w/phase shift
    CX = Y .* C2;
    GPSSignal = 2 * filtfilt(b, a, CX);
    GPSSignal = GPSSignal(padSize+1 : end-padSize);
end