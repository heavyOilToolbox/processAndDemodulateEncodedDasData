function binWords = decodeIrigBSignalToBits4(GPSSignal, t)
    
    % 3 sets of findpeaks parameters
    % first find marker bits (8 ms, 100 ms apart)
    % then look for bit markers (2 & 5 ms)

%     mpp = 0.0675;
    mpp = rms(GPSSignal);
    % parameters for data markers
    mpp2 = mpp;
    minpw2 = 0.006;
    maxpw2 = 0.011;
    mpd2 = 0.005;
    
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

    % first find data markers
    [~,loc2,~,~] = findpeaks(GPSSignal, t,...
        'MinPeakProminence', mpp2,...
        'MinPeakWidth', minpw2,...
        'MaxPeakWidth', maxpw2,...
        'MinPeakDistance', mpd2);
    
    % find zero bits
    [~,loc0,~,~] = findpeaks(GPSSignal, t,...
        'MinPeakProminence', mpp0,...
        'MinPeakWidth', minpw0,...
        'MaxPeakWidth', maxpw0,...
        'MinPeakDistance', mpd0);
    
    % find one bits
    [~,loc1,~,~] = findpeaks(GPSSignal, t,...
        'MinPeakProminence', mpp1,...
        'MinPeakWidth', minpw1,...
        'MaxPeakWidth', maxpw1,...
        'MinPeakDistance', mpd1);
    
%     if (size(loc2, 1) > size(loc2, 2))
%         loc = [loc2; loc0; loc1];
%     else
%         loc = [loc2, loc0, loc1]';
%     end
    loc = [loc2(:); loc0(:); loc1(:)];
    
    val = [2*ones(1,length(loc2)),...
           zeros(1,length(loc0)),...
           ones(1,length(loc1))];
    [~, sortIdx] = sort(loc);
    val = val(sortIdx);
    binWords = sprintf('%i', val);

end