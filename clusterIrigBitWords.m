function [IDX, d, meanSumD] = clusterIrigBitWords(irigBitWords)
warning off 
sigmaScale = 0.5;
N = length(irigBitWords);
C = nchoosek(1:N,2); %combinations for affinity matrix diagonal

%affinity matrix part 2: compute square euclidian distance
dX = zeros(size(C,1),1);
for kPair = 1:length(C)
    s = irigBitWords{C(kPair,1)};
    t = irigBitWords{C(kPair,2)};
    dxLev = LevenDistance(s, t) .^ 2;
    dxLen = (length(s) - 100) .^ 2 + (length(t) - 100) .^ 2;
    dX(kPair) = dxLev + dxLen;
end

    function L = gaussKernel(dX, C, sigmaScale)
        %affinity matrix part 3: gaussian kernel
        k = exp(-dX/2/(sigmaScale)^2); %gaussian kernel
        A = sparse(C(:,1), C(:,2), k, N, N);
        A = A + A';
        
        %construct diagonal matrix D whos diagonal are the equal to the row sums of
        %A (affinity matrix)
        D = sum(A, 2);
        
        %construct matrix L = D^(-1/2)*A*D^(-1/2)
        D = sqrt(D); D = 1./D; D = diag(D);
        L = D*A*D;
        L(isnan(L)) = 0;
    end

    function Y = kLargestEigenVectors(L, K)
        K = round(K);
        [X,D] = eig(L,'vector');
        [~,sortIdx] = sort(D);
        X = X(:,sortIdx);
        Y = X(:,(end-K+1):end);
        Y = bsxfun(@rdivide,Y,sum(abs(Y),2)); %normalize the rows of Y to length 1
        Y(~isfinite(Y)) = 0;
    end

    function IDX = spectralCluster(Y)
        K = size(Y, 2);
        %find the K most orthogonal eigenvectors
        dotProductY = abs(Y*Y.');
        dotProductY(logical(eye(length(Y)))) = nan;
        centroidIdx = randi(length(Y));
        dotProductY(centroidIdx,:) = nan; %dont pick the same row twice
        dotProductY = abs(dotProductY); %use the absolute value because looking for closer to zero
        for kCentroid = 2:K
            [~,im] = min(sum(dotProductY(:,centroidIdx),2)); %find next best centroid
            centroidIdx(kCentroid) = im; %reassign
            dotProductY(centroidIdx,:) = nan; %dont pick the same row twice
        end
        %cluster eigenvectors with kmeans, initialized with these centroids
        IDX = kmeans(abs(Y),[],'Start',abs(Y(centroidIdx,:)),'EmptyAction','drop','OnlinePhase','on');
        IDX(isnan(IDX)) = 0;
    end

        L = gaussKernel(dX, C, sigmaScale);
        d = sparse(C(:,1), C(:,2), dX, N, N);
        d = d + d';
        d = full(d);
        sumd = sum(d, 2);
        numCluster = 2;
        interClusterDistanceSum = inf;
        for kCluster = 2:6
            Y = kLargestEigenVectors(L, kCluster);
            idx = spectralCluster(Y);
            thisSumD = 0;
            for kc = 1:kCluster
                kDistance = sumd(idx==kc);
                meanDistance = mean(kDistance);
                intraClusterSqDistance = sum((kDistance - meanDistance) .^ 2);
%                 fprintf('kCluster %0.0f, kc %0.0f, dist %0.1f\n', kCluster, kc, intraClusterSqDistance);
                thisSumD = thisSumD + intraClusterSqDistance;
            end
%             fprintf('kCluster %0.0f, sumd %0.1f \n', kCluster, thisSumD);
            if thisSumD < interClusterDistanceSum
                numCluster = kCluster;
                interClusterDistanceSum = thisSumD;
                IDX = idx;
            end
        end
%         fprintf('numCluster %0.0f\n', numCluster);
        
        meanSumD = zeros(numCluster, 1);
        sumd = sum(d, 2);
        for kCluster = 1:numCluster
            meanSumD(kCluster) = mean(sumd(IDX==kCluster));
%             fprintf('IDX==%0.0f mean(sumd) %0.1f\n', kCluster, meanSumD(kCluster));
        end
        

%         sumD = sum(full(d), 2);

end