%% offDiagToCol takes in a variable amount of matrices, finds their upper off diagonal elements, 
%   and combines all respective elements and their corresponding cell pairs into a column vector
function col = offDiagToCol(varargin)
    % n represents col length
    n = 0;
    % prev keeps track of prev index for column vector concatenation
    prev = 0;
    % calculate length of column vector for preallocation
    for i=1:length(varargin)
        m = size(varargin{i},2);
        n = n + width(varargin{i})*length(varargin{i}) - m*(m+1)/2;
    end
    % col is initialized to empty matrix of n rows and three columns 
    % col1 contains values, col2 and col3 contain cell pair corresponding to value in col1
    col = zeros([n,3]);
    % select elements of each matrix's upper off diagonal and append to column vector 
    for j=1:length(varargin)
        temp = getIdxUpTri(varargin{j});
        % append to column vector
        col(prev+1:prev+height(temp),:) = temp;
        % update prev
        prev = prev + height(temp);
    end 
end