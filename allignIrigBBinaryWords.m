function allignedArray = allignIrigBBinaryWords(binaryWord, doVerbose)
% expected input binaryWord 1x28 cell of binary words
% expected output allignedArray 28 x 90 array of 1 and 0
    
    allignedArray = zeros(length(binaryWord), 90);
    for kWord = 1:length(binaryWord)
        
        irigBitWords = strsplit(binaryWord{kWord},'2');
        emptyMask = cellfun(@isempty,irigBitWords);
        [irigBitWords(emptyMask)] = [];
        nWord = length(irigBitWords);

        for k = 1:nWord
            irigBitWords{k} = pad(irigBitWords{k}, 9, 'right', '0');
            if length(irigBitWords{k}) > 9
                irigBitWords{k} = irigBitWords{k}(1:9);
            end
        end
        
        % force straight binary seconds to sit at the end if there are some
        % missing '2' separators
        if nWord < 10
            sbs1 = irigBitWords{end-1};
            sbs2 = irigBitWords{end};
            irigBitWords{9} = sbs1;
            irigBitWords{10} = sbs2;
            emptyMask = cellfun(@isempty, irigBitWords);
            if any(emptyMask)
                [irigBitWords{emptyMask}] = deal('000000000');
            end
        end
        if nWord > 10
            irigBitWords = [irigBitWords(1:8), irigBitWords(end-1:end)];
        end
        
        if doVerbose
            fprintf('%s ', irigBitWords{:});
            fprintf('\n');
        end
        
        thisString = pad(cell2mat(irigBitWords), 90, 'right', '0');
        if length(thisString) > 90
            fprintf('Irregular Binary Word %s\n', binaryWord{kWord});
            thisString = thisString(1:90);
        end
        allignedArray(kWord,:) = thisString - '0';
    end
end