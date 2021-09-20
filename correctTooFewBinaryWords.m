function binaryWord = correctTooFewBinaryWords(binaryWord)
    
    s = strsplit(binaryWord, '2');
    emptyMask = cellfun(@isempty, s);
    [s(emptyMask)] = [];
    wordLength = cellfun(@length, s);
    [m, im] = max(wordLength);
    if m < 14
        return;
    end
    if (im == 1)
        s = [s{1}(1:8), s{1}(9:end), s(2:end)];
    else
        if (m == 18)
            if (im == 9)
                s = [s(1:im-1) , s{im}(1:9), s{im}(10:18)];
            else
                s = [s(1:im-1) , s{im}(1:9), s{im}(10:18), s(im+1:end)];
            end
        elseif (m == 19)
            if (im == 9)
                s = [s(1:im-1) , s{im}(1:9), s{im}(11:19)];
            else
                s = [s(1:im-1) , s{im}(1:9), s{im}(11:19), s(im+1:end)];
            end
        else
            n1 = floor(m / 2);
            if (im == 9)
                s = [s(1:im-1) , s{im}(1:n1), s{im}(n1+1:end)];
            else
                s = [s(1:im-1) , s{im}(1:n1), s{im}(n1+1:end), s(im+1:end)];
            end
        end
    end
    binaryWord = correctTooFewBinaryWords(['2', strjoin(s, '2'), '2']);
end
