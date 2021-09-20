function [s, m, h, d, y, tos, bs] = decodeIrigBSignalBits(irigBitWords)
    P = zeros(10, 9);
    P(1, :) = [1,   2,    4,    8,    0,    10,    20,    40,    0];   % seconds
    P(2, :) = [1,   2,    4,    8,    0,    10,    20,    40,    0];   % minutes. 9th bit unused
    P(3, :) = [1,   2,    4,    8,    0,    10,    20,    0,     0];   % hours.   8th & 9th bits unused
    P(4, :) = [1,   2,    4,    8,    0,    10,    20,    40,    80];  % day.
    P(5, :) = [100, 200,  0,    0,    0,    0.1,   0.2,   0.4,   0.8]; % day (cont) & tenths of seconds
    P(6, :) = [1,   2,    4,    8,    0,    10,    20,    40,    80];  % year (00-99)
    P(7, :) = [0,   0,    0,    0,    0,    0,     0,     0,     0];   % TODO: device control function
    P(8, :) = [0,   0,    0,    0,    0,    0,     0,     0,     0];   % TODO: device control function
    P(9, :) = [1,   2,    4,    8,    16,   32,    64,    128,   256]; % straight binary seconds (0-86399)
    P(10,:) = [512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 0];   % straight binary seconds (cont)
    
    irigBitWords = strsplit(irigBitWords,'2');
    emptyMask = cellfun(@isempty,irigBitWords);
    [irigBitWords(emptyMask)] = [];
    
    nWord = length(irigBitWords);
    if nWord > 0
        secondBits = irigBitWords{1} - '0';
    else
        secondBits = zeros(1,8);
    end
    if nWord > 1
        minuteBits = irigBitWords{2} - '0';
    else
        minuteBits = zeros(1,9);
    end
    if nWord > 2
        hourBits = irigBitWords{3} - '0';
    else
        hourBits = zeros(1,9);
    end
    if nWord > 3
        dayBits1 = irigBitWords{4} - '0';
    else
        dayBits1 = zeros(1,9);
    end
    if nWord > 4
        dayBits2 = irigBitWords{5} - '0';
    else
        dayBits2 = zeros(1,9);
    end
    if nWord > 8
        binarySecondBits1 = irigBitWords{9} - '0';
    else
        binarySecondBits1 = zeros(1,9);
    end
    if nWord > 9
        binarySecondBits2 = irigBitWords{10} - '0';
    else
        binarySecondBits2 = zeros(1,9);
    end
    
    nb = length(secondBits);
    if nb < 8
        secondBits = [secondBits, zeros(1, 8-nb)];
    end
    nb = length(minuteBits);
    if nb < 9
        minuteBits = [minuteBits, zeros(1, 9-nb)];
    end
    nb = length(hourBits);
    if nb < 7
        hourBits = [hourBits, zeros(1, 7-nb)];
    end
    nb = length(dayBits1);
    if nb < 9
        dayBits1 = [dayBits1, zeros(1, 9-nb)];
    end
    nb = length(dayBits2);
    if nb < 2
        dayBits2 = [dayBits2, zeros(1, 2-nb)];
    end
    nb = length(binarySecondBits1);
    if nb < 9
        binarySecondBits1 = [binarySecondBits1, zeros(1, 9-nb)];
    end
    nb = length(binarySecondBits2);
    if nb < 9
        binarySecondBits2 = [binarySecondBits2, zeros(1, 9-nb)];
    end
    
    s = sum(secondBits(1:8) .* P(1, 1:8));
    m = sum(minuteBits(1:9) .* P(2, 1:9));
    h = sum(hourBits(1:7) .* P(3, 1:7));
    d = sum([dayBits1(1:9) .* P(4, 1:9), dayBits2(1:2) .* P(5, 1:2)]);
    y = 0;
    tos = 0;
    bs = sum([binarySecondBits1(1:9) .* P(9, :), binarySecondBits2(1:8) .* P(10, 1:8)]);
   
    if (s > 60)
        s = -1;
    end
    if (m > 60)
        m = -1;
    end
    if (h > 24)
        h = -1;
    end
    if (d > 365)
        d = -1;
    end
    if (bs > 86400)
        bs = -1;
    end
    
end