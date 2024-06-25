%% function 'buildRatMatrix' takes in cell array C and builds n by 8 matrix 'all', 
%  where n is number of elements in upper right off-diagonal matrix for
%  each statistical calculation (cov, r-vals, etc.) totaled together from
%  each recording. The columns of 'all' consist of the following:
% col1: cov, col2: r-vals, col3: p-vals, col4: L95CI, col5: U95CI, col6:
% recording number, col7: cell 1, col8: cell2.

function all = buildRatMatrix(C)
    n = 0;
    % calculate length of column vector for preallocation
    for i=1:length(C)
         m = size(C{i},2);
         n = n + width(C{i})*length(C{i}) - m*(m+1)/2;
    end
    all = zeros([n,8]);
    % prev keeps track of prev index for column vector concatenation
    prev = 0;
    % loop traverses through each recording
    % length(C) captures num of recordings
    for i = 1:length(C)
        % find cov, rVals, pVals, L95CI, U95CI of each cell recording
        % store upper right off-diagonal of data into column vectors
        aCov = offDiagToCol(cov([C{i,:}]));
        [rVals,pVals,L95CI,U95CI] = corrcoef([C{i,:}]);
        rVals = offDiagToCol(rVals);
        pVals = offDiagToCol(pVals);
        L95CI = offDiagToCol(L95CI);
        U95CI = offDiagToCol(U95CI);
        % store recording number
        recNum = repmat(i,height(aCov),1);
        % matrix of data for recording i 
        temp = [aCov(:,1) rVals(:,1) pVals(:,1) L95CI(:,1) U95CI(:,1) recNum aCov(:,2) aCov(:,3)];
        % append recording data to 'all' 
        all(prev+1:prev+height(temp),:) = temp;
        % update prev
        prev = prev + height(temp);
    end
end