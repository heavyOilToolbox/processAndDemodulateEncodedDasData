function GPSSignal = andreasEnvelopeDetectorIrigB(IRIGBSignal, fc, fs)
    % andreas' envelope detector
    cut_off = (fc / 2) / (fs / 2); % Cut-off frequency: ~500Hz
    [b, a] = butter(1, cut_off, 'low');
    t = 0:1/fs:(length(IRIGBSignal) / fs - 1 / fs); t = t';
    I = IRIGBSignal .* sin(2 * pi * t * fc);
    Q = IRIGBSignal .* cos(2 * pi * t * fc);
    I_filt = filtfilt(b, a, I);
    Q_filt = filtfilt(b, a, Q);
    GPSSignal = 2 * sqrt(I_filt.*I_filt + Q_filt.*Q_filt);
end