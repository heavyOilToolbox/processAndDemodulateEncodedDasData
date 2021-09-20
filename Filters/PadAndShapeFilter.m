function Y = PadAndShapeFilter(X,fs)
    % make pulse shaper
    filterRolloff = 0.21;
    filterSpanPulses = 10;
    samplesPerPeriod = 10;
    b = rcosdesign(filterRolloff, filterSpanPulses, samplesPerPeriod, 'normal');
    b = b / sum(b);
    padSize = sqrt(fs);
    X = padarray(X, [padSize, 0], 0, 'both');
    Y = filtfilt(b, 1, X);
    Y = Y(padSize+1:end-padSize);
end