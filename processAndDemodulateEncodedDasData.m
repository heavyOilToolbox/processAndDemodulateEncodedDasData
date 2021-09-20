function processAndDemodulateEncodedDasData
% TODO: use mcc instead of mcr to compile
% mcc -m -I /IrigB -d /IrigB/processAndDemodulateEncodedDasData processAndDemodulateEncodedDasData
% TODO: parameters....
% profile on;
dasDataDirectory = input('Enter DAS dPhase Root Directory ', 's');
% pathParts = strsplit(dasDataDirectory, filesep);
% rootDriveLetter = pathParts{1};
if isdeployed
    startMarkerA = input('Enter Data Marker A Start Channel ', 's');
    stopMarkerA = input('Enter Data Marker A End Channel ', 's');
    startMarkerB = input('Enter Data Marker B Start Channel ', 's');
    stopMarkerB = input('Enter Data Marker B End Channel ', 's');
    startDataMarkerA = str2double(startMarkerA);
    stopDataMarkerA = str2double(stopMarkerA);
    startDataMarkerB = str2double(startMarkerB);
    stopDataMarkerB = str2double(stopMarkerB);
else
    startDataMarkerA = 93;
    stopDataMarkerA = 97;
    startDataMarkerB = 145;
    stopDataMarkerB = 150;
end

% delta time vector
fc = 1000; % Hz. IRIG-B carrier frequency
% nPulse = 16384; % TODO: read from binary
% dt = (0:(nPulse - 1)) / fs; % seconds
scaleFactor = 10430; 

% list binary files
if ~exist('dd','var')
    dd = dirRecursive3(dasDataDirectory, '\.bin$');
end
nFile = length(dd);
fprintf('%0.0f Files Found\n', nFile);

% [header, zones, depthCal, nominalDepth, measuredDepth, frameData, qualityData] = readPinnacleDasBinary(d{1});
[fileHeader, zones, ~, ~, ~, ~, ~] = readPinnacleDasBinary(dd{1});

fs = fileHeader.SamplingRate;
%nPulse = double(header.NumberOfFrames);
nPulse = double(fileHeader.FrameCapacity);
nLocation = double(fileHeader.NumberOfChannels);
quality_block_size = double(fileHeader.QualityBlockSize);
channelList = [];
nZone = zones(1);
for kZone = 1:nZone
    startChannel = zones(2 + 3 * (kZone - 1));
    endChannel = zones(3 + 3 * (kZone - 1));
    stride = zones(4 + 3 * (kZone - 1));
    channelList = [channelList, startChannel:stride:endChannel];
end
% first file timestamp / start of acquisition
acquisitionStartTime = datetime(fileHeader.Year,...
                            fileHeader.Month,...
                            fileHeader.Day,...
                            fileHeader.Hour,...
                            fileHeader.Minute,...
                            fileHeader.Second,...
                            double(fileHeader.Microsecond) / 1000);

startA = find(channelList==startDataMarkerA);
stopA = find(channelList==stopDataMarkerA);
startB = find(channelList==startDataMarkerB);
stopB = find(channelList==stopDataMarkerB);
if isempty(startA)
    error('Data Marker A Start Channel %i Out of Range', startDataMarkerA);
end
if isempty(stopA)
    error('Data Marker A End Channel %i Out of Range', stopDataMarkerA);
end
if isempty(startB)
    error('Data Marker B Start Channel %i Out of Range', startDataMarkerB);
end
if isempty(stopB)
    error('Data Marker B End Channel %i Out of Range', stopDataMarkerB);
end

% data marker A / fiber stretcher #1: IRIG-B 120 timecode from symmetricom
piezo1mask = false(nLocation, 1);
piezo1mask(startA:stopA) = true;

% data marker B / fiber stretcher #2: also IRIG-B 120 timecode from symmetricom
piezo2mask = false(nLocation, 1);
piezo2mask(startB:stopB) = true;

piezoMask = piezo1mask | piezo2mask;

% data tables for excel output
T1 = table({}, [], [], [], [], [], [], [], [], [], [], [], [], [], {}, ...
    'VariableNames', {'FileName',...
    'Year', 'Month', 'Day', 'Hour', 'Minute', 'Second', 'Microsecond', 'BinarySecond',... 
    'InterpIRIGBinarySecond', 'DeltaTimeMicrosecond', 'StartPulse', 'PulseCount', 'ElapsedPulsesClock', 'FullFilePath'});

T2 = table([], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], ...
    'VariableNames', {'SecondsMarkerPulse',...
        'BestGuessSecond', 'BestGuessMinute', 'BestGuessHour', 'BestGuessBinarySecond',...
        'IRIGBSecond', 'IRIGBCalcSecond',...
        'IRIGBMinute', 'IRIGBCalcMinute',...
        'IRIGBHour', 'IRIGBCalcHour',...
        'Year', 'IRIGBCalcMonth', 'IRIGBCalcDay', ...
        'IRIGBBinarySeconds', 'IRIGBCalcBinarySeconds', 'IRIGBDayOfYear'});

nFilePerChunk = lcm(nPulse, fs) / nPulse; % e.g. 625
nDataChunk = floor(nFile / nFilePerChunk) + 1;
nPulsePreviousChunk = 0;
for kChunk = 1:nDataChunk
% for kChunk = nDataChunk:nDataChunk % debug tail
    fileTimeStamp = zeros(nFilePerChunk, 7);
    fileDateTime = NaT(nFilePerChunk, 1);
    elapsedPulsesClock = zeros(nFilePerChunk, 1);
    filePulseCount = zeros(nFilePerChunk, 1);
    %fprintf('Processing Data Chunk %3.0f\n', kChunk);
    piezoIrigWeightedAvg = zeros(nPulse * nFilePerChunk, 1);
    lastValidRow = nFilePerChunk * nPulse;
    for kFile = 1:nFilePerChunk
        fileIdx = kFile + (kChunk - 1) * nFilePerChunk;
        firstRow = 1 + (kFile - 1) * nPulse;
        lastRow = nPulse * kFile;
        if (fileIdx == nFile)
            lastValidRow = lastRow;
        end
        if (fileIdx > nFile)
            piezoIrigWeightedAvg(firstRow:lastRow) = 0;
            continue;
        end
        binaryFileName = dd{fileIdx};
        fprintf('Reading File %03i/%03i (Chunk %i) %06i/%06i (Total) %s\n', kFile, nFilePerChunk, kChunk, fileIdx, nFile, binaryFileName);
        % [~, x, quality] = readPinnacleDasBinary(binaryFileName);
        [fileHeader, dPhase, quality] = readPinnacleDasBinary(binaryFileName, piezoMask);
        
        % process timestamps
        fileTimeStamp(kFile,:) = [double(fileHeader.Year),...
            double(fileHeader.Month),...
            double(fileHeader.Day),...
            double(fileHeader.Hour),...
            double(fileHeader.Minute),...
            double(fileHeader.Second),...
            double(fileHeader.Microsecond)];
        
        thisFileDateTime = datetime(fileHeader.Year,...
            fileHeader.Month,...
            fileHeader.Day,...
            fileHeader.Hour,...
            fileHeader.Minute,...
            fileHeader.Second,...
            double(fileHeader.Microsecond) / 1000);
        fileDateTime(kFile) = thisFileDateTime;
        
        % stack the stretcher channels
        piezoIrigWeightedAvg(firstRow:lastRow, 1) = stackSignal(dPhase, quality, nPulse, quality_block_size);
        
        % estimate number of pulses ? (may not always be 16384)
        elapsedPulsesClock(kFile) = fs * seconds(thisFileDateTime - acquisitionStartTime);
        filePulseCount(kFile) = double(fileHeader.NumberOfFrames);
        
    end
    
    w = cumsum(piezoIrigWeightedAvg) / scaleFactor;
    firstNanIdx = length(w);
    hasNaN = false;
    if (lastValidRow < lastRow)
        firstNanIdx = find(isnan(w), 1, 'first');
        if ~isempty(firstNanIdx)
            hasNaN = true;
            w(firstNanIdx:end) = 0;
        else
            w(lastValidRow+1:end) = 0;
        end
    end
    
    fprintf('Detrending Stacked IRIG-B Signal\n');
    % use highpass filter to detrend signal
    B = fir1(nPulse, 10/5000, 'high');
    if (lastValidRow < lastRow)
        if (hasNaN)
            wf = fftDetrendSignal(B, w(1:firstNanIdx-1));
            %wf = fastIIRFilter(B, w(1:firstNanIdx-1));
        else
            wf = fftDetrendSignal(B, w(1:lastValidRow));
            %wf = fastIIRFilter(B, w(1:lastValidRow));
        end
    else
        wf = fftDetrendSignal(B, w);
        %wf = fastIIRFilter(B, w);
    end
    
    % fprintf('Demodulating AM Signal from Chunk %0.0f\n', kChunk);
    % [PD1, PD2, PD3, PD4, PD5, PD6,...
    %     PD1f, PD2f, PD3f, PD4f, PD5f, PD6f,...
    %     fullMedianSignal, fullMedianSignalShaped] = StackProductDetectors(w, fc, fs, initialPhase, length(w), false);
    %[PD1, PD2, PD3, PD4, PD5, fullMedianSignal] = StackProductDetectorsBasic(w, fc, fs, initialPhase, length(w), false);
    
    fprintf('Demodulating AM Signal from Chunk %0.0f (Multiply and Shift Detector)\n', kChunk);
    if isdeployed
        PD2 = MultiplyAndShiftDetector(wf, fc, fs);
    else
        PD2 = gather(gpuMultiplyAndShiftDetector(wf, fc, fs));
    end
    fprintf('Demodulating AM Signal from Chunk %0.0f (High Bandwidth Lock In Detector)\n', kChunk);
    if isdeployed
        PD4 = HighBWLockInDetector(wf, fc, fs);
    else
        PD4 = gather(gpuHighBWLockInDetector(wf, fc, fs));
    end
    fprintf('Locating Seconds Markers in IRIG-B Timecode (xcorrFindSecondsMarker)\n');
    %allLags = gather(gpuXcorrFindSecondsMarker(wf, fc, fs, nPulse, false));
    allLags = xcorrFindSecondsMarker(gather(wf), fc, fs, nPulse, false);
    if (allLags(1) > 0)
        secondsMarkerPulse = [allLags(1)-fs; allLags];
    else
        secondsMarkerPulse = allLags;
    end
    dataFileName = sprintf('%s\\DataChunk%06i.mat', dasDataDirectory, kChunk);
    
    % split demodulated IRIG-B signal at lag positions from allLags
    nLag = length(allLags);
    binaryWords = cell(nLag + 1, 6);
    PD = PD2 + PD4;
    [S, M, H, D, BS, CS, CM, CH, CBS] = deal(zeros(nLag+1,1));
    for k = 0:nLag
        if (k == 0)
            startIdx = 1;
            endIdx = min(allLags(1), fs);
            nLeftPad = fs - endIdx;
            z1 = cat(1, zeros(nLeftPad, 1), PD(startIdx:endIdx));
            z2 = cat(1, zeros(nLeftPad, 1), PD2(startIdx:endIdx));
            z3 = cat(1, zeros(nLeftPad, 1), PD4(startIdx:endIdx));
        elseif (k == nLag)
            startIdx = allLags(end);
            endIdx = length(PD);
            nRightPad = fs - (endIdx-startIdx+1);
            if (nRightPad > 0)
                z1 = cat(1, PD(startIdx:endIdx), zeros(nRightPad, 1));
                z2 = cat(1, PD2(startIdx:endIdx), zeros(nRightPad, 1));
                z3 = cat(1, PD4(startIdx:endIdx), zeros(nRightPad, 1));
            else
                endIdx = startIdx + fs - 1;
                z1 = PD(startIdx:endIdx);
                z2 = PD2(startIdx:endIdx);
                z3 = PD4(startIdx:endIdx);
            end
        else
            startIdx = allLags(k);
            endIdx = startIdx + fs - 1;
            z1 = PD(startIdx:endIdx);
            z2 = PD2(startIdx:endIdx);
            z3 = PD4(startIdx:endIdx);
        end
        %fprintf('Decoding Irig-B Signal from %s to bits\n', dataFileName);
        binaryWords{k+1, 1} = decodeIrigBSignalToBits(z1);
        binaryWords{k+1, 2} = decodeIrigBSignalToBits(z2);
        binaryWords{k+1, 3} = decodeIrigBSignalToBits(z3);
        dt = (0:(fs+20)-1) / fs; dt = dt';
        binaryWords{k+1, 4} = decodeIrigBSignalToBits4(padarray(z1, [10,0], 0, 'both'), dt);
        binaryWords{k+1, 5} = decodeIrigBSignalToBits4(padarray(z2, [10,0], 0, 'both'), dt);
        binaryWords{k+1, 6} = decodeIrigBSignalToBits4(padarray(z3, [10,0], 0, 'both'), dt);
        fprintf('%s\n',binaryWords{k+1,:})
        %fprintf('Alligning Binary Words\n');
        try
            allignedString = allignIrigBBinaryWords(binaryWords(k+1,:), false);
            [s, m, h, d, ~, ~, bs] = decodeIrigBSignalBitsFromArray(mode(allignedString,1));
        catch
            s = 0; m = 0; h = 0; d = 0; bs = 0;
        end
        calculatedSecond = mod(bs, 60);
        calculatedMinute = mod(bs - calculatedSecond, 3600) / 60;
        calculatedHour = (bs - calculatedSecond - calculatedMinute * 60) / 3600;
        calculatedBinarySeconds = s + m * 60 + h * 3600;
        doVerbose = true;
        if doVerbose
            fprintf('%3.0f %3.0f %3.0f %4.0f %6.0f | ', s, m, h, d, bs);
            fprintf('%3.0f %3.0f %3.0f     %6.0f\n', calculatedSecond, calculatedMinute, calculatedHour, calculatedBinarySeconds);
        end
        S(k+1) = s;
        M(k+1) = m;
        H(k+1) = h;
        D(k+1) = d;
        BS(k+1) = bs;
        CS(k+1) = calculatedSecond;
        CM(k+1) = calculatedMinute;
        CH(k+1) = calculatedHour;
        CBS(k+1) = calculatedBinarySeconds;
    end
    
    currentYear = mode(fileTimeStamp(:,1));
    dv = datevec(datenum(double(currentYear), 1, D, 0, 0, 0));
    YMD = dv(:,1:3);
    
    startFileIdx = 1 + (kChunk - 1) * nFilePerChunk;
    endFileIdx = min(nFile, kChunk * nFilePerChunk);
    nTableRow = endFileIdx - startFileIdx + 1;
    thisFileList = dd(startFileIdx:endFileIdx);
    [~, dasFName, ~] = cellfun(@fileparts, thisFileList, 'uni', 0);
    
    %filePulseCount = 1 + cumsum(filePulseCount) - filePulseCount(1);
    filePulseCount = filePulseCount(1:length(thisFileList));
    
    % time correction
    deltaBS = [diff(BS); -1];
    validDeltaBS = any(deltaBS == [1, -86399], 2);
    deltaCBS = [diff(CBS); -1];
    validDeltaCBS = any(deltaCBS == [1, -86399], 2);
    deltaS = [diff(S); -1];
    validDeltaS = any(deltaS == [1, -59], 2);
    deltaCS = [diff(CS); -1];
    validDeltaCS = any(deltaCS == [1, -59], 2);
    deltaM = [diff(M); -1];
    validDeltaM = any(deltaM == [0, 1, -59], 2);
    deltaCM = [diff(CM); -1];
    validDeltaCM = any(deltaCM == [0, 1, -59], 2);
    deltaH = [diff(H); -1];
    validDeltaH = any(deltaH == [0, 1, -23], 2);
    deltaCH = [diff(CH); -1];
    validDeltaCH = any(deltaCH == [0, 1, -23], 2);
    
    validBinarySecondMask = (BS == CBS) & validDeltaBS;
    validCalcBinarySecondMask = (BS == CBS) & validDeltaCBS;
    
    validSecondMask = (S == CS) & ( validDeltaS & (deltaS == deltaCS) );
    validCalcSecondMask = validBinarySecondMask & validDeltaCS;
    
    validMinuteMask = (M == CM) & validDeltaM;
    validCalcMinuteMask = validBinarySecondMask & validDeltaCM;
    
    validHourMask = (H == CH) & validDeltaH;
    validCalcHourMask = validBinarySecondMask & validDeltaCH;

    corrS = zeros(length(secondsMarkerPulse), 1);
    m1 = validSecondMask & validCalcSecondMask;
    m2 = ~validSecondMask & validCalcSecondMask;
    m3 = validSecondMask & ~validCalcSecondMask;
    %m4 = ~(m1 | m2 | m3);
    validCorrSecondsMask = (m1 | m2 | m3);
    corrS(m1) = S(m1);
    corrS(m2) = CS(m2);
    corrS(m3) = S(m3);
    % interpolation fails where seconds wrap around!
    %interpS = interp1(secondsMarkerPulse(~m4), corrS(~m4), secondsMarkerPulse(m4), 'linear', 'extrap');
    
    corrM = zeros(length(secondsMarkerPulse), 1);
    m1 = validMinuteMask & validCalcMinuteMask;
    m2 = ~validMinuteMask & validCalcMinuteMask;
    m3 = validMinuteMask & ~validCalcMinuteMask;
    %m4 = ~(m1 | m2 | m3);
    validCorrMinutesMask = (m1 | m2 | m3);
    corrM(m1) = M(m1);
    corrM(m2) = CM(m2);
    corrM(m3) = M(m3);
    
    corrH = zeros(length(secondsMarkerPulse), 1);
    m1 = validHourMask & validCalcHourMask;
    m2 = ~validHourMask & validCalcHourMask;
    m3 = validHourMask & ~validCalcHourMask;
    %m4 = ~(m1 | m2 | m3);
    validCorrHoursMask = (m1 | m2 | m3);
    corrH(m1) = H(m1);
    corrH(m2) = CH(m2);
    corrH(m3) = H(m3);
    
    corrBS = zeros(length(secondsMarkerPulse), 1);
    m1 = validBinarySecondMask & validCalcBinarySecondMask;
    m2 = ~validBinarySecondMask & validCalcBinarySecondMask;
    m3 = validBinarySecondMask & ~validCalcBinarySecondMask;
    validCorrBSMask = (m1 | m2 | m3);
    corrBS(m1) = BS(m1);
    corrBS(m2) = CBS(m2);
    corrBS(m3) = BS(m3);
    
    % backfill gaps in corrected seconds with corrBS where valid
    backFillSecondsMask = ~validCorrSecondsMask & validCorrBSMask;
    backFillMinutesMask = ~validCorrMinutesMask & validCorrBSMask;
    backFillHoursMask = ~validCorrHoursMask & validCorrBSMask;
    if any(backFillSecondsMask)
        backFillSeconds = mod(corrBS(backFillSecondsMask), 60);
        corrS(backFillSecondsMask) = backFillSeconds;
        validCorrSecondsMask = validCorrSecondsMask | backFillSecondsMask;
    end
    if any(backFillMinutesMask)
        backFillMinutes = mod(corrBS(backFillMinutesMask) - mod(corrBS(backFillMinutesMask), 60), 3600) / 60;
        corrM(backFillMinutesMask) = backFillMinutes;
        validCorrMinutesMask = validCorrMinutesMask | backFillMinutesMask;
    end
    if any(backFillHoursMask)
        backFillHours = floor(corrBS(backFillHoursMask) / 3600);
        corrH(backFillHoursMask) = backFillHours;
        validCorrHoursMask = validCorrHoursMask | backFillHoursMask;
    end
    
    % backfill gaps in binary seconds with a calculated value where valid
    backFillBSMask = validCorrSecondsMask & validCorrMinutesMask & validCorrHoursMask & ~validCorrBSMask;
    if any(backFillBSMask)
        backFillBS = corrS(backFillBSMask) + corrM(backFillBSMask) * 60 + corrH(backFillBSMask) * 3600;
        corrBS(backFillBSMask) = backFillBS;
        validCorrBSMask = validCorrBSMask | backFillBSMask;
    end
    
    %interpolate remaining gaps in binary seconds
    if any(validCorrBSMask) && any(~validCorrBSMask)
        interpBS = interp1(secondsMarkerPulse(validCorrBSMask), corrBS(validCorrBSMask), secondsMarkerPulse(~validCorrBSMask), 'linear', 'extrap');
        corrBS(~validCorrBSMask) = round(interpBS);
    end
    
    [X, iu] = unique(secondsMarkerPulse);
    if (length(iu) > 1)
        V = corrBS(iu);
    else
        V = [];
    end
    
    % FINALLY ... interpolate remaining gaps in seconds, minutes, hours
    if any(~validCorrSecondsMask) && any(validCorrSecondsMask) && ~isempty(V)
        interpBS = round(interp1(X, V, secondsMarkerPulse(~validCorrSecondsMask), 'linear', 'extrap'));
        interpS = mod(interpBS, 60);
        corrS(~validCorrSecondsMask) = interpS;
    end
    if any(~validCorrMinutesMask) && any(validCorrMinutesMask) && ~isempty(V)
        interpBS = round(interp1(X, V, secondsMarkerPulse(~validCorrMinutesMask), 'linear', 'extrap'));
        interpM = floor(mod(interpBS, 3600) / 60);
        corrM(~validCorrMinutesMask) = interpM;
    end
    if any(~validCorrHoursMask) && any(validCorrHoursMask) && ~isempty(V)
        interpBS = round(interp1(X, V, secondsMarkerPulse(~validCorrHoursMask), 'linear', 'extrap'));
        interpH = floor(interpBS / 3600);
        corrH(~validCorrHoursMask) = interpH;
    end

    %corrBS = corrS + corrM * 60 + corrH * 3600;
    
    fileTimeStamp = fileTimeStamp(1:nTableRow, :);
    fileHour = fileTimeStamp(:, 4);
    fileMinute = fileTimeStamp(:, 5);
    fileSecond = fileTimeStamp(:, 6);
    fileMicro = fileTimeStamp(:, 7);
    fileBinarySecond = fileHour * 3600 + fileMinute * 60 + fileSecond + fileMicro / 1000000;
    
    % estimate correct timestamp by interpolation
    chunkPulseCount = cumsum([1;filePulseCount]); chunkPulseCount = chunkPulseCount(1:end-1);
    
    [xs, iu] = unique(secondsMarkerPulse);
    interpBinarySecond = interp1(xs, corrBS(iu), chunkPulseCount, 'linear', 'extrap');
    deltaBinarySecond = 1000000 * (interpBinarySecond - fileBinarySecond);
    
    % offset pulse counts by sum of previous chunk
    chunkPulseCount = chunkPulseCount + nPulsePreviousChunk;
    % secondsMarkerPulse = secondsMarkerPulse + (kChunk - 1) * nPulse * nFilePerChunk;
    secondsMarkerPulse = secondsMarkerPulse + nPulsePreviousChunk;
    nPulsePreviousChunk = nPulsePreviousChunk + sum(filePulseCount);
    
    t1 = table(dasFName', fileTimeStamp(1:nTableRow,1), fileTimeStamp(1:nTableRow,2),...
        fileTimeStamp(1:nTableRow,3), fileHour, fileMinute, ...
        fileSecond, fileMicro, fileBinarySecond,...
        interpBinarySecond, deltaBinarySecond, chunkPulseCount, filePulseCount, elapsedPulsesClock(1:nTableRow), thisFileList', ...
        'VariableNames', {'FileName',...
        'Year', 'Month', 'Day', 'Hour', 'Minute', 'Second', 'Microsecond', 'BinarySecond',... 
        'InterpIRIGBinarySecond', 'DeltaTimeMicrosecond', 'StartPulse', 'PulseCount', 'ElapsedPulsesClock', 'FullFilePath'});
    
    t2 = table(secondsMarkerPulse, corrS, corrM, corrH, corrBS,...
        S, CS, M, CM, H, CH,...
        YMD(:,1), YMD(:,2), YMD(:,3),...
        BS, CBS, D, ...
        'VariableNames', {'SecondsMarkerPulse',...
        'BestGuessSecond', 'BestGuessMinute', 'BestGuessHour', 'BestGuessBinarySecond',...
        'IRIGBSecond', 'IRIGBCalcSecond',...
        'IRIGBMinute', 'IRIGBCalcMinute',...
        'IRIGBHour', 'IRIGBCalcHour',...
        'Year', 'IRIGBCalcMonth', 'IRIGBCalcDay', ...
        'IRIGBBinarySeconds', 'IRIGBCalcBinarySeconds', 'IRIGBDayOfYear'});
    
    % append data table
    T1 = [T1; t1];
    T2 = [T2; t2];
    
    fprintf('Saving %s\n', dataFileName);
    save(dataFileName, 'w', 'PD2', 'PD4', 'allLags', 'binaryWords', 't1', 't2');
    
%     % debug only
%     xlFileName = sprintf('foo%i.xlsx', kChunk);
%     writetable(t1, xlFileName, 'Sheet', 1);
%     writetable(t2, xlFileName, 'Sheet', 2);
    
%     
%     kNugget = 0;
%     startIdx = 1 + kNugget * fs;
%     endIdx = fs;
%    
%     binaryWords = cell(nFilePerChunk * nPulse / fs, 28);
% 
%     kWord = 0;
% 
%     while (endIdx <= length(w))
%         kWord = kWord + 1;
% 
%         medianSignal = fullMedianSignal(startIdx:endIdx);
% %         medianSignalShaped = fullMedianSignalShaped(startIdx:endIdx);
% 
%         binaryWords{kWord,1} = decodeIrigBSignalToBits(medianSignal);
% %         binaryWords{kWord,2} = decodeIrigBSignalToBits(medianSignalShaped);
%         binaryWords{kWord,3} = decodeIrigBSignalToBits(PD1(startIdx:endIdx));
%         binaryWords{kWord,4} = decodeIrigBSignalToBits(PD2(startIdx:endIdx));
%         binaryWords{kWord,5} = decodeIrigBSignalToBits(PD3(startIdx:endIdx));
%         binaryWords{kWord,6} = decodeIrigBSignalToBits(PD4(startIdx:endIdx));
%         binaryWords{kWord,7} = decodeIrigBSignalToBits(PD5(startIdx:endIdx));
% %         binaryWords{kWord,8} = decodeIrigBSignalToBits(PD6(startIdx:endIdx));
% %         binaryWords{kWord,9} = decodeIrigBSignalToBits(PD1f(startIdx:endIdx));
% %         binaryWords{kWord,10} = decodeIrigBSignalToBits(PD2f(startIdx:endIdx));
% %         binaryWords{kWord,11} = decodeIrigBSignalToBits(PD3f(startIdx:endIdx));
% %         binaryWords{kWord,12} = decodeIrigBSignalToBits(PD4f(startIdx:endIdx));
% %         binaryWords{kWord,13} = decodeIrigBSignalToBits(PD5f(startIdx:endIdx));
% %         binaryWords{kWord,14} = decodeIrigBSignalToBits(PD6f(startIdx:endIdx));
%         
%         X1 = padarray(medianSignal, [10,0], 0, 'both');
% %         X2 = padarray(medianSignalShaped, [10,0], 0, 'both');
%         X3 = padarray(PD1(startIdx:endIdx), [10,0], 0, 'both');
%         X4 = padarray(PD2(startIdx:endIdx), [10,0], 0, 'both');
%         X5 = padarray(PD3(startIdx:endIdx), [10,0], 0, 'both');
%         X6 = padarray(PD4(startIdx:endIdx), [10,0], 0, 'both');
%         X7 = padarray(PD5(startIdx:endIdx), [10,0], 0, 'both');
% %         X8 = padarray(PD6(startIdx:endIdx), [10,0], 0, 'both');
% %         X9 = padarray(PD1f(startIdx:endIdx), [10,0], 0, 'both');
% %         X10 = padarray(PD2f(startIdx:endIdx), [10,0], 0, 'both');
% %         X11 = padarray(PD3f(startIdx:endIdx), [10,0], 0, 'both');
% %         X12 = padarray(PD4f(startIdx:endIdx), [10,0], 0, 'both');
% %         X13 = padarray(PD5f(startIdx:endIdx), [10,0], 0, 'both');
% %         X14 = padarray(PD6f(startIdx:endIdx), [10,0], 0, 'both');
%         % pad signal with zeros
%         t2 = (0:length(X1)-1) / fs; t2 = t2';
%         
%         binaryWords{kWord,15} = decodeIrigBSignalToBits4(X1, t2);
%         binaryWords{kWord,16} = decodeIrigBSignalToBits4(X2, t2);
%         binaryWords{kWord,17} = decodeIrigBSignalToBits4(X3, t2);
%         binaryWords{kWord,18} = decodeIrigBSignalToBits4(X4, t2);
%         binaryWords{kWord,19} = decodeIrigBSignalToBits4(X5, t2);
%         binaryWords{kWord,20} = decodeIrigBSignalToBits4(X6, t2);
%         binaryWords{kWord,21} = decodeIrigBSignalToBits4(X7, t2);
%         binaryWords{kWord,22} = decodeIrigBSignalToBits4(X8, t2);
%         binaryWords{kWord,23} = decodeIrigBSignalToBits4(X9, t2);
%         binaryWords{kWord,24} = decodeIrigBSignalToBits4(X10, t2);
%         binaryWords{kWord,25} = decodeIrigBSignalToBits4(X11, t2);
%         binaryWords{kWord,26} = decodeIrigBSignalToBits4(X12, t2);
%         binaryWords{kWord,27} = decodeIrigBSignalToBits4(X13, t2);
%         binaryWords{kWord,28} = decodeIrigBSignalToBits4(X14, t2);
% 
%         fprintf('%i\n',kWord);
%         fprintf('%s\n',binaryWords{kWord,:})
% 
%         kNugget = kNugget + 1;
%         startIdx = 1 + kNugget * fs;
%         endIdx = startIdx + fs - 1;
%     end
%     
%     binaryWordFileName = sprintf('%s\\DataChunk%06iBinaryWords.mat', rootDriveLetter, kChunk);
%     fprintf('Saving %s\n', binaryWordFileName);
%     save(binaryWordFileName, 'binaryWords')
    
end
    % save spreadsheet and data table to disk
    pathParts = strsplit(dasDataDirectory, filesep);
    sessionFolderName = pathParts{end-1};
    xlFileName = sprintf('%s.xlsx', sessionFolderName);
    xlFullFileName = fullfile(dasDataDirectory, xlFileName);
    writetable(T1, xlFullFileName, 'Sheet', 1);
    writetable(T2, xlFullFileName, 'Sheet', 2);
    matFileName = sprintf('%s.mat', sessionFolderName);
    matFullFileName = fullfile(dasDataDirectory, matFileName);
    save(matFullFileName, 'T1', 'T2');
    
% profile viewer; %debug
end

function y = stackSignal(x, quality, nPulse, qualityBlockSize)
    % % make array for weighted average based on std and avg quality
    % % 16064 x 1 -> 64 x 251
    avg_quality = quality(1:4:end);
    % std_quality = quality(4:4:end);
    avgQuality = reshape(avg_quality, nPulse / qualityBlockSize, []); % avg ./ std
    %avgQuality = avgQuality(:, piezoMask);
    avgQuality = avgQuality(:);
    avgQuality = repmat(avgQuality', qualityBlockSize, 1);
    nCol = size(x, 2);
    avgQuality = reshape(avgQuality(:), nPulse, nCol);
    avgQuality = avgQuality ./ sum(avgQuality, 2);
    y = sum(double(x) .* avgQuality, 2); %weighted average avg/std
end

% this hangs on while executing on deployed GPUs



%     tks = regexp(dasFName, '(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})(\d{2})_(\d+)', 'tokens');
%     tks = [tks{:}];
%     ntk = length(tks);
%     fdv = zeros(ntk, 7);
%     for k = 1:ntk
%         tk = tks{k};
%         fdv(k,1) = str2double(tk{1});
%         fdv(k,2) = str2double(tk{2});
%         fdv(k,3) = str2double(tk{3});
%         fdv(k,4) = str2double(tk{4});
%         fdv(k,5) = str2double(tk{5});
%         fdv(k,6) = str2double(tk{6});
%         fdv(k,7) = str2double(tk{7});
%     end