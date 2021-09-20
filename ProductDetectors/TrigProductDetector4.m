function GPSSignal = TrigProductDetector4(IRIGBSignal, fc, fs, varargin)
% trig identity: sin^2(t) * cos^2(t) = (1 = cos(4*t)) / 8
%     Z = TrigProductDetector4(Y,Fc,Fs) demodulates the amplitude modulated signal Y from
%     the carrier frequency Fc (Hz). Y and Fc have sample frequency Fs (Hz).
%     The modulated signal Y has zero initial phase, and zero carrier
%     amplitude, for suppressed carrier modulation. A lowpass filter is used
%     in the demodulation.  The default filter is: [NUM,DEN] =
%     butter(5,Fc*2/Fs).
%     
%     Z = TrigProductDetector4(Y,Fc,Fs,INI_PHASE) specifies the initial phase (rad) of the
%     modulated signal.
%  
%     Z = TrigProductDetector4(Y,Fc,Fs,INI_PHASE,CARRAMP) specifies the carrier amplitude
%     of the modulated signal for transmitted carrier modulation.
%  
%     Z = TrigProductDetector4(Y,Fc,Fs,INI_PHASE,CARRAMP,NUM,DEN) specifies the filter to
%     be used in the demodulation. 

    [IRIGBSignal, fc, fs, INI_PHASE, ~, b, a] = ParseInputArguments(IRIGBSignal, fc, fs, varargin{:});
    
    initialPhase = pi * INI_PHASE;

    padSize = sqrt(fs);
    [dt, Y] = detrendAndPadSignal(IRIGBSignal, padSize, fs);
    C2 = sin(2 * pi * fc * dt + initialPhase); % carrier w/phase shift
    C3 = cos(2 * pi * fc * dt + initialPhase); % out-of-phase carrier
    W = Y .* C2 .* C3 .* C3; % product of signal and carrier and oop carrier^2
    GPSSignal = 8 * filtfilt(b, a, W); % lowpass filter to remove carrier
    GPSSignal = GPSSignal(padSize+1 : end-padSize);
end