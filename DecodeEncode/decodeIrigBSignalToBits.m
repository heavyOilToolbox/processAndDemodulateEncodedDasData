function binWords = decodeIrigBSignalToBits(GPSSignal)
    
    % 3 sets of findpeaks parameters
    % first find marker bits (8 ms, 100 ms apart)
    % then look for bit markers (2 & 5 ms)
    mpp = rms(GPSSignal);
%     mpp = 0.0675;

    % parameters for data markers 8 ms
    mpp2 = mpp;
    minpw2 = 0.006;
    maxpw2 = 0.011;
    mpd2 = 0.06;
    
    % parameters for `zero` bits
    mpp0 = mpp;
    minpw0 = 0.001;
    maxpw0 = 0.0035;
    mpd0 = 0.005;
    
    % parameters for `one` bits
    mpp1 = mpp;
    minpw1 = maxpw0;
    maxpw1 = minpw2;
    mpd1 = 0.004;
    
    % TODO: eliminate magic numbers
    bw = cell(10,1);
    for k = 1:10
        startIdx = 1 + (k - 1) * 1000;
        endIdx = k * 1000;
%         X = GPSSignal(startIdx:endIdx);
%         fprintf('startIdx: %i, endIdx: %i, N: %i\n', startIdx, endIdx, length(GPSSignal));
        X = padarray(GPSSignal(startIdx:endIdx), [10, 0], 0, 'both');
        dt = -0.001:0.0001:0.1009;
        % find zero bits
        [~,loc0k,~,~] = findpeaks(X, dt,...
            'MinPeakProminence', mpp0,...
            'MinPeakWidth', minpw0,...
            'MaxPeakWidth', maxpw0,...
            'MinPeakDistance', mpd0);
        % find one bits
        [~,loc1k,~,~] = findpeaks(X, dt,...
            'MinPeakProminence', mpp1,...
            'MinPeakWidth', minpw1,...
            'MaxPeakWidth', maxpw1,...
            'MinPeakDistance', mpd1);
        if (k == 1)
            loc2k = [0; 0.1];
        else
            loc2k = [0.1];
        end
        loc = [loc2k(:); loc0k(:); loc1k(:)];
    
        val = [2*ones(1,length(loc2k)),...
               zeros(1,length(loc0k)),...
               ones(1,length(loc1k))];
        [~, sortIdx] = sort(loc);
        val = val(sortIdx);
        bw{k} = sprintf('%i', val);
    end
    binWords = [bw{:}];

end