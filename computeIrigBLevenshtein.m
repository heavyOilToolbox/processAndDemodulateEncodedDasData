function L = computeIrigBLevenshtein(irigBitWords)
% irigBitWords: cell array of IRIG B strings
% L: lower triangular matrix of levenshteitn string distances
    nString = length(irigBitWords);
    L = zeros(nString);
    for k1 = 1:nString
        s1 = irigBitWords{k1};
        for k2 = 1:nString
            if k1 >= k2
                continue;
            end
            s2 = irigBitWords{k2};
            L(k1, k2) = lev(s1, s2);
            L(k2, k1) = L(k1, k2);
        end
    end
end