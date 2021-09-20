function [S, M, H, D, BS, CS, CM, CH, CBS] = allignAndDecodeIrigBString(binaryWords, doVerbose)
    nChunk = size(binaryWords, 1);
    [S, M, H, D, BS, CS, CM, CH, CBS] = deal(zeros(nChunk,1));
    for kChunk = 1:nChunk
        bw = binaryWords(kChunk, :);
        [IDX, ~, meanSumD] = clusterIrigBitWords(bw);
        %IDX(1) is the cluster containing the expected string
        clusterMask = (IDX == IDX(1));
        
        if sum(clusterMask) == 1
            [~, sortIdx] = sort(meanSumD);
            cluster1Mask = IDX==sortIdx(1);
            cluster2Mask = IDX==sortIdx(2);
            clusterMask = cluster1Mask | cluster2Mask;
        end
        clusterMask(1) = 0; %exclude expected timestamp duh
        bw = bw(clusterMask);
        
        % correct for too many/too few binary words
        for k = 1:length(bw)
            bw{k} = correctTooFewBinaryWords(bw{k});
        end
        
        allignedString = allignIrigBBinaryWords(bw, false);
        [s, m, h, d, ~, ~, bs] = decodeIrigBSignalBitsFromArray(mode(allignedString,1));
        
        calculatedSecond = mod(bs, 60);
        calculatedMinute = mod(bs - calculatedSecond, 3600) / 60;
        calculatedHour = (bs - calculatedSecond - calculatedMinute * 60) / 3600;
        calculatedBinarySeconds = s + m * 60 + h * 3600;
        
        if doVerbose
            fprintf('%3.0f %3.0f %3.0f %4.0f %6.0f | ', s, m, h, d, bs);
            fprintf('%3.0f %3.0f %3.0f     %6.0f\n', calculatedSecond, calculatedMinute, calculatedHour, calculatedBinarySeconds);
        end
       S(kChunk) = s;
       M(kChunk) = m;
       H(kChunk) = h;
       D(kChunk) = d;
       BS(kChunk) = bs;
       CS(kChunk) = calculatedSecond;
       CM(kChunk) = calculatedMinute;
       CH(kChunk) = calculatedHour;
       CBS(kChunk) = calculatedBinarySeconds;
    end
end