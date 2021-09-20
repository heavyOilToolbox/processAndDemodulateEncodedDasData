function Y = fastIIRFilter(b, X)
    phaseShiftPoints = length(b);
    halfPhaseShiftPoints = floor(phaseShiftPoints / 2);
% %     X0 = gpuArray(padarray(X,[halfPhaseShiftPoints,0],'both'));
% %     X0 = gpuArray(padarray(X,[halfPhaseShiftPoints,0], 'symmetric', 'both'));
%     nPoint = length(X);
%     midPoint = floor(nPoint / 2);
%     X1 = gpuArray(X(1:midPoint + halfPhaseShiftPoints));
%     X2 = gpuArray(X((midPoint - halfPhaseShiftPoints + 1):end));
% %     X1 = X0(1:(midPoint + halfPhaseShiftPoints));
% %     X2 = X0((midPoint - halfPhaseShiftPoints + 1):nPoint);
%     Y1 = fftfilt(b, X1); % forward
%     Y2 = flipud(fftfilt(b, flipud(X2)));
%     Y = gather([Y1(halfPhaseShiftPoints+1:end); Y2(1:end-halfPhaseShiftPoints)]);

    % detrend the head and tail segments then pad X with the same trend
    x = [(1:phaseShiftPoints)', ones(phaseShiftPoints,1)];
    startPad = X(1:phaseShiftPoints);
    if feature('IsDebugMode')
        fprintf('fastIIRFilter: Computing Linear Trend on Left Half of Signal... ');
    end
    p1 = x \ startPad;
    postPad = X(end-phaseShiftPoints+1:end);
    if feature('IsDebugMode')
        fprintf('Right Half of Signal\n');
    end
    p2 = x \ postPad;
    
    if feature('IsDebugMode')
        fprintf('fastIIRFilter: Computing Detrended Left Pad Values... ');
    end
    startPadDetrended = startPad - x * p1;
    if feature('IsDebugMode')
        fprintf('Right Pad Values\n');
    end
    endPadDetrended = postPad - x * p2;
    
    if feature('IsDebugMode')
        fprintf('fastIIRFilter: Computing Left Pad x Values... ');
    end
    xPrePad = x - [phaseShiftPoints, 0];
    if feature('IsDebugMode')
        fprintf('Right Pad x Values\n');
    end
    xPostPad = x + [phaseShiftPoints, 0];
    %xPostPad = [((halfPhaseShiftPoints+1):phaseShiftPoints-1)', ones(halfPhaseShiftPoints,1)];
    if feature('IsDebugMode')
        fprintf('fastIIRFilter: Computing Left Pad y Values... ');
    end
    yPrePad =  startPadDetrended + xPrePad * p1;
    if feature('IsDebugMode')
        fprintf('Right Pad y Values\n');
    end
    yPostPad = endPadDetrended + xPostPad * p2;
    
    % debug graphics
    %subplot(2,2,1)
    %plot(x(:,1), startPad, x(:,1), x * p1, xPrePad, yPrePad);
    %subplot(2,2,2)
    %plot(x(:,1), startPadDetrended);
    %subplot(2,2,3)
    %plot(x(:,1), postPad, x(:,1), x * p2, xPostPad, yPostPad);
    %subplot(2,2,4)
    %plot(x(:,1), endPadDetrended);
    
    if feature('IsDebugMode')
        fprintf('fastIIRFilter: Concatenating [yPrePad; X; yPostPad]\n');
    end
    cpuX0 = [yPrePad; X; yPostPad];
    if feature('IsDebugMode')
        fprintf('fastIIRFilter: Moving Padded Array to GPU\n');
    end
    X0 = gpuArray(cpuX0);
    %X0 = gpuArray([yPrePad; X; yPostPad]);

    nPoint = length(X0);
    midPoint = floor(nPoint / 2);

    
    % X1 = X0(1:(midPoint + phaseShiftPoints + halfPhaseShiftPoints)); % nope?
    leftEndIdx = midPoint + halfPhaseShiftPoints;
    if feature('IsDebugMode')
        fprintf('fastIIRFilter: Slicing gpuArray Left Side X1 = X0(%i:%i)... ', 1, leftEndIdx);
    end
    X1 = X0(1:leftEndIdx); % bingo
    % X2 = X0((midPoint - phaseShiftPoints - halfPhaseShiftPoints + 1):nPoint); % ?
    % X2 = X0((midPoint - phaseShiftPoints - halfPhaseShiftPoints + 1):end); %nope
    rightStartIdx = midPoint - halfPhaseShiftPoints + 1;
    rightEndIdx = length(X0);
    if feature('IsDebugMode')
        fprintf('Right Side X2 = X0(%i:%i)\n', rightStartIdx, rightEndIdx);
    end
    X2 = X0(rightStartIdx : rightEndIdx); % bingo
    
    if feature('IsDebugMode')
        fprintf('fastIIRFilter: Forward filtering pass Y1 = fftfilt(b, X1)... ');
    end
    Y1 = fftfilt(b, X1); % forward
    if feature('IsDebugMode')
        fprintf('Reverse filtering pass Y2 = flipud(fftfilt(b, flipud(X2)))\n');
    end
    Y2 = flipud(fftfilt(b, flipud(X2)));
    %Y = gather([Y1(phaseShiftPoints+1:end); Y2(1:end-phaseShiftPoints)]);
%     Y = [Y1(phaseShiftPoints+halfPhaseShiftPoints+1:end); Y2(1:end-phaseShiftPoints-halfPhaseShiftPoints)]; %nope
%     Y = [Y1(phaseShiftPoints+halfPhaseShiftPoints+1:end); Y2(halfPhaseShiftPoints+1:end-phaseShiftPoints)]; %nope
%     Y = [Y1(phaseShiftPoints+halfPhaseShiftPoints+1:end); Y2(phaseShiftPoints+1:end-halfPhaseShiftPoints)]; %nope
    y1StartIdx = phaseShiftPoints + halfPhaseShiftPoints + 1;
    y1EndIdx = length(Y1);
    y2StartIdx = 1;
    y2EndIdx = length(Y2) - phaseShiftPoints - halfPhaseShiftPoints;
    if feature('IsDebugMode')
        fprintf('fastIIRFilter: Concatenating Forward and Reverse Arrays: [Y1(%i:%i); Y2(%i:%i)]\n', y1StartIdx, y1EndIdx, y2StartIdx, y2EndIdx);
    end
    Y = [Y1(y1StartIdx:y1EndIdx); Y2(y2StartIdx:y2EndIdx)]; 
end