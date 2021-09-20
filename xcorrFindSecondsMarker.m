function allLags = xcorrFindSecondsMarker(wf, fc, fs, nPulsePerBlock, doPlot)
    % fs = 10000;
    % fc = 1000;
    % nPulse = 16384;
    nPulse = nPulsePerBlock;

    % make reference IRIG-B seconds marker
    secondMarkerPulseWidth = floor(8 * fs / fc);
    firstFrameSecondMarkerIdx = 1;
    firstFrameLastSecondMarkerStartIdx = round(990 * fs / fc) + 1;
    firstFrameLastSecondMarkerEndIdx = firstFrameLastSecondMarkerStartIdx + secondMarkerPulseWidth - 1;
    secondFrameSecondMarkerStartIdx = fs + 1;
    secondFrameSecondMarkerEndIdx = secondFrameSecondMarkerStartIdx + secondMarkerPulseWidth - 1;
    initialPhase = pi / 2;
    x = zeros(nPulse, 1);
    x(firstFrameSecondMarkerIdx:secondMarkerPulseWidth) = 1;
    x(firstFrameLastSecondMarkerStartIdx:firstFrameLastSecondMarkerEndIdx) = 1;
    x(secondFrameSecondMarkerStartIdx:secondFrameSecondMarkerEndIdx) = 1;
    y = ammod(x/6, fc, fs, initialPhase);
    %y = gpuArray(y);
    
    N = length(wf);
    amplitudeMax = max(abs(wf));
    nColZ = floor(N / secondFrameSecondMarkerEndIdx) - 1;

columnLag = zeros(nColZ, 1);
for k = 1:nColZ
    startIdx = 1 + (k-1) * secondFrameSecondMarkerEndIdx;
    endIdx = startIdx + 2 * secondFrameSecondMarkerEndIdx;
    zf = wf(startIdx:endIdx);
    amplitudeThreshold = rms(zf);
    
    zf(zf > amplitudeThreshold) = amplitudeMax;
    zf(zf < -amplitudeThreshold) = -amplitudeMax;
    
    % validDataStartIdx = nRowZ - 10080 + 1;
    % validDataEndIdx = 2 * nRowZ - 9901;
    % nValidPoint = validDataEndIdx - validDataStartIdx + 1;
    [C0, LAGS0] = xcorr(zf, y);
    % C0 = C0(validDataStartIdx:validDataEndIdx);
    % LAGS0 = LAGS0(validDataStartIdx:validDataEndIdx);
    
    [~, im] = max(abs(C0));
    columnLag(k) = LAGS0(im);
end
allLags = (secondFrameSecondMarkerEndIdx * (0:nColZ-1))' + columnLag;
if (doPlot)
    plot(wf);hold on;plot(allLags,wf(allLags),'o','MarkerFaceColor','r');hold off
end
    
    
%     %validDataStartIdx = length(wf) - 10080 + 1;
%     %validDataEndIdx = 2 * length(wf) - 9901;
%     
%     % find the first data marker
%     % amplify the peaks in the first set of nPulse by thresholding
%     %zf = reshape(zf, nPulse, []);
%     nTrailingFramePulse = nPulse - mod(nSample, nPulse);
%     zf = reshape(cat(1, wf, zeros(nTrailingFramePulse, 1, 'gpuArray')), nPulse, []);
%     amplitudeThreshold = rms(zf);
%     amplitudeMax = max(abs(wf));
%     zf(zf > amplitudeThreshold) = amplitudeMax;
%     zf(zf < -amplitudeThreshold) = -amplitudeMax;
%     
%     %zf = wf(1:nPulse);
%     %amplitudeThreshold = rms(zf);
%     %amplitudeMax = max(abs(zf));
%     %zf(zf > amplitudeThreshold) = amplitudeMax;
%     %zf(zf < -amplitudeThreshold) = -amplitudeMax;
%     [C0, LAGS0] = xcorr(zf(:,1), y);
%     [~, im] = max(abs(C0(1:nPulse-1)));
%     initialLag0 = LAGS0(im);
%     [~, im] = max(abs(C0(nPulse:end)));
%     initialLag1 = LAGS0(nPulse + im - 1);
%     if (abs(initialLag0 - initialLag1) > 1000)
%         initialLag = initialLag0;
%     else
%         initialLag = round((initialLag0 + (initialLag1 - fs)) / 2);
%     end
%     
%     nTrailingFramePulse = nPulse - mod(nSample, nPulse);
%     nColZ = size(zf, 2);
%     nRowY = 2 * nPulse;
%     nColY = floor(nColZ / 2);
%     yf = reshape(zf(:,1:2*nColY), nRowY, nColY);
%     columnLag = zeros(nColY, 2);
% %     columnLag = zeros(nColY, 5);
% %     columnLag2 = zeros(nColY, 1);
% %     columnLag3 = zeros(nColY, 1);
% %     columnLag4 = zeros(nColY, 1);
% %     columnLag5 = zeros(nColY, 1);
%     for k = 1:nColY
%         validDataStartIdx = nRowY - 10080 + 1;
%         validDataEndIdx = 2 * nRowY - 9901;
%         nValidPoint = validDataEndIdx - validDataStartIdx + 1;
%         [C0, LAGS0] = xcorr(yf(:,k), y);
%         C0 = C0(validDataStartIdx:validDataEndIdx);
%         LAGS0 = LAGS0(validDataStartIdx:validDataEndIdx);
%         [~, im] = max(abs(C0));
%         firstLag = LAGS0(im) + 2 * nPulse * (k - 1);
%         %columnLag(k, 3) = LAGS0(im) + 2 * nPulse * (k - 1);
%         
%         midPoint = floor(nValidPoint / 2);
%         if (im > midPoint)
%             [~, im1] = max(abs(C0(1:im-128)));
%             columnLag(k, 1) = LAGS0(im1) + 2 * nPulse * (k - 1);
%             columnLag(k, 2) = firstLag;
%         else
%             [~, im2] = max(abs(C0(im+127:end)));
%             columnLag(k, 1) = firstLag;
%             columnLag(k, 2) = LAGS0(im2 + im + 127 - 1) + 2 * nPulse * (k - 1);
%         end
%         
%         
% %         if ((im - 2 * fs) > 0)
% %             startIdx = max(1, (im - 2 * fs - 1024));
% %             endIdx = min(2*nPulse, im - 2 * fs + 1023);
% %             [m1, im1] = max(abs(C0(startIdx:endIdx)));
% %             columnLag(k, 1) = LAGS0(im1 + startIdx - 1) + 2 * nPulse * (k - 1);
% %         end
% %         if ((im - fs) > 0)
% %             startIdx = max(1, (im - fs - 1024));
% %             endIdx = min(2*nPulse, im - fs + 1023);
% %             [m2, im2] = max(abs(C0(startIdx:endIdx)));
% %             columnLag(k, 2) = LAGS0(im2 + startIdx - 1) + 2 * nPulse * (k - 1);
% %         end
% %         if ((im + fs) < nValidPoint)
% %             startIdx = max(1, (im + fs - 1023));
% %             endIdx = min(2 * nPulse, im + fs + 1024);
% %             [m3, im3] = max(abs(C0(startIdx:endIdx)));
% %             columnLag(k, 4) = LAGS0(im3 + startIdx - 1) + 2 * nPulse * (k - 1);
% %         end
% %         if ((im + 2 * fs) < nValidPoint)
% %             startIdx = max(1, (im + 2 * fs - 1023));
% %             endIdx = min(2 * nPulse, im + 2 * fs + 1024);
% %             [m4, im4] = max(abs(C0(startIdx:endIdx)));
% %             columnLag(k, 5) = LAGS0(im4 + startIdx - 1) + 2 * nPulse * (k - 1);
% %         end
%         
%             
%         
%         
%     end
%     allColumnLag = ravel(columnLag');
%     allColumnLag = allColumnLag(allColumnLag~=0);
%     lagSecondsMarkerIdx = floor(allColumnLag / fs);
%     badIncrementMask = (diff(lagSecondsMarkerIdx)==0);
%     allColumnLag(badIncrementMask) = [];
%     lagSecondsMarkerIdx(badIncrementMask) = [];
%     
%     allSecondsMarkerIdx = (0:(length(wf) / fs))';
%     % lagSecondsMarkerIdx = floor(allColumnLag / fs);
%     missingSecondsMarkerIdx = ~ismember(allSecondsMarkerIdx, lagSecondsMarkerIdx);
%     interpSecondsMarker = interp1(lagSecondsMarkerIdx, allColumnLag, allSecondsMarkerIdx(missingSecondsMarkerIdx), 'linear', 'extrap');
%     
%     allLags = zeros(size(allSecondsMarkerIdx));
%     allLags(missingSecondsMarkerIdx) = round(interpSecondsMarker);
%     allLags(~missingSecondsMarkerIdx) = allColumnLag;
%     allLags(allLags > length(wf)) = [];
%     
%     if (doPlot)
%         figure; plot(wf); hold on; plot(allLags, wf(allLags), 'o', 'MarkerFaceColor', 'r'); hold off;
%     end
%     %plot(interpSecondsMarker, wf(interpSecondsMarker), 'o', 'MarkerFaceColor', 'y');
    
%     
%     % plot(LAGS0, abs(C0), [initialLag0,initialLag1], [m0, m1], 'o', 'MarkerFaceColor', 'r')
%     
%     % cross correlate the signals
% %     zf = zf(1:end-nTrailingFramePulse);
% %     [C, LAGS] = xcorr(zf, y);
%     [C, LAGS] = xcorr(zf(1:end-nTrailingFramePulse), y);
%     %idx0 = nSample - nPulse + 1;
%     %C = abs(C);
%     C = C(validDataStartIdx:validDataEndIdx);
%     LAGS = LAGS(validDataStartIdx:validDataEndIdx);
%     
%     % there are 10080 lags less than zero
%     %firstBlockC = C(LAGS <= 0);
%     %firstBlockLag = LAGS(LAGS <= 0);
%     %[initialLagVal, im] = max(abs(C(LAGS <= 0)));
%     %initialLag = LAGS(im);
%     
%     %[initialLagVal, im] = max(C(idx0:nSample));
%     %initialLag = LAGS(idx0 + im - 1);
%     
% 
%     %allLags = zeros(nDataFrame, 1, 'gpuArray');
%     %allLagVal = zeros(nDataFrame, 1, 'gpuArray');
%     %allLags(1) = initialLag;
%     %allLagVal(1) = initialLagVal;
% 
%     %LAGS = LAGS(nSample:end);
%     %C = C(nSample:end);
%     
%     % there are length(w) + 180 - 1 LAGS greater than zero
%     %C = abs(C(LAGS > 0));
%     %LAGS = LAGS(LAGS > 0);
% 
%     timeWindowSamples = 2048;
% %     timeWindowSamples = 8192;
%     startIdx = max(1, fs + initialLag - timeWindowSamples / 2);
%     endIdx = min(fs, startIdx + timeWindowSamples - 1);
%     actualTimeWindowSamples = endIdx - startIdx + 1;
%     
%     % reshape C for parallellization
%     %nCorrDataFrame = floor(length(C) / fs) - 1;
%     nCorr = length(C);
%     nCorrDataFrame = ceil(nCorr / fs);
%     nTrailingCorrSample = fs - mod(nCorr, fs);
%     
%     timeWindowMask = false(fs, nCorrDataFrame, 'gpuArray');
%     timeWindowMask(startIdx:endIdx, :) = true;
%     B = abs(cat(2, C, zeros(1, nTrailingCorrSample, 'gpuArray')));
%     [m, im] = max(reshape(B(timeWindowMask), actualTimeWindowSamples, nCorrDataFrame));
%     allStartIdx = fs * (1:nCorrDataFrame) + initialLag - timeWindowSamples / 2;
%     %allLags1 = fs + [initialLag, LAGS(im + allStartIdx) - 1]';
%     %allLagVal1 = [initialLagVal, m]';
%     lagIdx = im + allStartIdx;
%     validLagMask = lagIdx <= nCorr;
%     lagIdx = lagIdx(validLagMask);
%     allLags1 = fs + (LAGS(lagIdx) - 1)';
%     allLagVal1 = m(validLagMask)';
%     %allLags1(1) = initialLag + fs;
%     
%     
%     % replaced with logical mask above
%     % B = reshape(C(1:fs*nCorrDataFrame), 10000, nCorrDataFrame);
%     % if nTrailingCorrSample
%     %     B(:, end+1) = zeros(fs, 1, 'gpuArray');
%     %     B(1:nTrailingCorrSample, end) = C(end-nTrailingCorrSample+1:end);
%     % end
%     % B = B(startIdx:endIdx, :);
% 
%     % deprecated by parallellized code above
%     %for k = 1:(nDataFrame - 1)
%     %    startIdx = fs * k + initialLag - timeWindowSamples / 2;
%     %    endIdx = startIdx + timeWindowSamples - 1;
%     %    [m, im] = max(C(startIdx:endIdx));
%     %    thisLag = LAGS(im + startIdx - 1);
%     %    % fprintf('{%0.0f,%0.0f} thisLag: %0.0f\n', startIdx, endIdx, thisLag);
%     %    allLags(k+1) = thisLag;
%     %    allLagVal(k+1) = m;
%     %end
% 
% %     % check for trailing seconds marker if the numerator of nDataFrame is not
% %     % integral
% %     if (nTrailingPulse)
% %         k = nDataFrame - 1; %-1?
% %         startIdx = max(1, fs * (k + 1) + initialLag - 1000);
% %         endIdx = min(length(C), startIdx + 2000 - 1);
% %         if (startIdx < nSample) && ((startIdx + 1080) < nSample)
% %             [m, im] = max(C(startIdx:endIdx));
% %             thisLag = LAGS(im + startIdx - 1);
% %             allLags1(end+1) = thisLag + fs;
% %             allLagVal1(end+1) = m;
% %         end
% %     end
% 
%     if doPlot
%         figure; plot(LAGS, C, allLags1-fs, allLagVal1, 'o')
%         xlabel('Lag (Pulses)'); ylabel('Cross-Correlation');
%         figure; plot(wf); hold on; plot(allLags1, wf(allLags1), 'o', 'MarkerFaceColor', 'r'); hold off
%         xlabel('Pulses'); ylabel('IRIG-B Signal (rad)');
%     end
% 
%     % allLags = allLags + fs;
%     allLags = allLags1;
end