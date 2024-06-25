%% function 'buildRatMatrixByObXPc' takes in two cell arrays and builds n by 8 matrix 'all', 
%  where n is number of elements in upper right rectangular matrix (ob x pc) for
%  each statistical calculation (cov, r-vals, etc.) totaled together from
%  each recording. The columns of 'all' consist of the following:
% col1: cov, col2: r-vals, col3: p-vals, col4: L95CI, col5: U95CI, col6:
% recording number, col7: cell 1, col8: cell2.

function all = buildRatMatrixByObXPc(varargin)
    n = 0;
    % calculate length of matrix 'all' for preallocation
    for i=1:length(varargin(:,:))
         m = size(varargin{i},2);
         n = n + width(varargin{i})*length(varargin{i}) - m*(m+1)/2;
    end
    all = zeros([n,8]);
    % prev keeps track of prev index for column vector concatenation
    prev = 0;
    % loop traverses through each recording
    % can use varargin{1} to capture num of recordings since ob and pc will
    %   have same num of recordings for each drug state
    for i = 1:length(varargin{1})
        % store nOb and nPc for calculations
        nOb = width(varargin{1}{i,:});
        nPc = width(varargin{2}{i,:});
        % find cov, rVals, pVals, L95CI, U95CI of ob x pc
        % store upper right rectangle (obXpc) of data into column vectors
        aCov = getRect(cov([varargin{1}{i,:} varargin{2}{i,:}]),nOb,nPc);
        [rVals,pVals,L95CI,U95CI] = corrcoef([varargin{1}{i,:} varargin{2}{i,:}]);
        rVals = getRect(rVals,nOb,nPc);
        pVals = getRect(pVals,nOb,nPc);
        L95CI = getRect(L95CI,nOb,nPc);
        U95CI = getRect(U95CI,nOb,nPc);
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