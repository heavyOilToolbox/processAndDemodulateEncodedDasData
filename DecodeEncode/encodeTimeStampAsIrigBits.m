function irigBString = encodeTimeStampAsIrigBits(s, m, h, d)

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
    
    S = encodeSeconds(s);
    M = encodeMinutes(m);
    H = encodeHours(h);
    [D1, D2] = encodeDays(d);
    bs = s + m * 60 + h * 3600;
    %fprintf('bs: %0.0f\n', bs);
    BS = encodeBinarySeconds(bs);
    
%     A = [
%         [S,0];...
%         M;...
%         H;...
%         D1;...
%         D2;...
%         zeros(2,9);...
%         BS(1:9);...
%         BS(10:18)
%         ];
    irigBString = [ '2',...
        sprintf('%i', S), '2',...
        sprintf('%i', M), '2',...
        sprintf('%i', H), '2',...
        sprintf('%i', D1), '2',...
        sprintf('%i', D2), '2',...
        sprintf('%i', zeros(1,9)), '2',...
        sprintf('%i', zeros(1,9)), '2',...
        sprintf('%i', zeros(1,9)), '2',...
        sprintf('%i', BS(1:9)), '2',...
        sprintf('%i', BS(10:18)), '2'];
end

function BS = encodeBinarySeconds(bs)
    BS = zeros(1,18);
    for k = 16:-1:1
        [q, bs] = divmod(bs, 2^k);
        if (q == 1)
            BS(k+1) = 1;
        end
    end
    if (bs == 1)
        BS(1) = 1;
    end
end

function [D1, D2] = encodeDays(d)
    D1 = zeros(1,9);
    D2 = zeros(1,9);
    [q, d] = divmod(d, 200);
    if (q == 1)
        D2(2) = 1;
    end
    [q, d] = divmod(d, 100);
    if (q == 1)
        D2(1) = 1;
    end
    [q, d] = divmod(d, 80);
    if (q == 1)
        D1(9) = 1;
    end
    D1(1:8) = encodeSeconds(d);
end

function H = encodeHours(h)
    H = zeros(1, 9);
    [q, h] = divmod(h, 20);
    if (q == 1)
        H(7) = 1;
    end
    [q, h] = divmod(h, 10);
    if (q == 1)
        H(6) = 1;
    end
    [q, h] = divmod(h, 8);
    if (q == 1)
        H(4) = 1;
    end
    [q, h] = divmod(h, 4);
    if (q == 1)
        H(3) = 1;
    end
    [q, h] = divmod(h, 2);
    if (q == 1)
        H(2) = 1;
    end
    if (h == 1)
        H(1) = h;
    end
end

function S = encodeSeconds(s)
    S = zeros(1, 8);
    [q, s] = divmod(s, 40);
    if (q == 1)
        S(8) = 1;
    end
    h = encodeHours(s);
    S(1:7) = h(1:7);
end

function M = encodeMinutes(m)
    M = [encodeSeconds(m), 0];
end