function [s, m, h, d, y, tos, bs] = decodeIrigBSignalBitsFromArray(irigBitWords)
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
    
    Q = ravel(P')';
    
    x = irigBitWords .* Q;
    s = sum(x(:, 1:8), 2);
    m = sum(x(:, 10:17), 2);
    h = sum(x(:, 19:26), 2);
    d = sum(x(:, 27:36), 2) + sum(x(:, 37:38), 2);
    y = sum(x(:, 46:54), 2);
    tos = zeros(size(x,1), 1);
    bs = sum(x(:, 73:81), 2) + sum(x(:, 82:90), 2);
    
    s(s > 60) = -1;
    m(m > 60) = -1;
    h(h > 24) = -1;
    d(d > 365) = -1;
    bs(bs > 86400) = -1;

end