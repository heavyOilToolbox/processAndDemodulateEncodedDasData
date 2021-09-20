function y = fftDetrendSignal(B, x)
    n = length(x);
    if (isdeployed || (n < 131072))
        [bb,aa] = butter(1, 10/5000, 'high');
        y = detrend(filtfilt(bb, aa, x));
        %y = gpuArray(z);
    else
        y = fastFIRFilter(B, x);
        t = gpuArray([(1:n)', ones(n,1)]);
        p = t(1:2:n,:) \ y(1:2:n); % decimate?
        y = y - t * p;
    end
end