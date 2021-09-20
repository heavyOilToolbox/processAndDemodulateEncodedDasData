function [t, Y] = detrendAndPadSignal(X, N, fs)
    Y = padarray(detrend(X), [N,0], 0, 'both');
    t = (0:length(Y)-1) / fs; t = t';
end