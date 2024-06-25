%% getIdxUpTri finds upper off-diagonal triangular matrix and returns 'v' containing values and their respective cell pairs
function v = getIdxUpTri(A)
    % store original matrix A as 'temp' for later calculations requiring original matrix size    
    temp = A;
    zeroIdxs = temp == 0;
    temp(zeroIdxs) = eps;
    % find linear indices of upper off-diagonal triangular matrix
    idx = ind2sub(size(temp),find(triu(temp,1)));
    % get original values that linear indicies 'idx' correspond to and store in 'A'
    A = A(idx);
    restoreZeroIdxs = A == eps;
    A(restoreZeroIdxs) = 0;
    % store corresponding cell pairs (c1 = 'r', c2 = 'c')
    [r,c] = ind2sub(size(temp),idx);
    % 'v' contains upper off-diagonal triangular matrix values as 'A', and
    %   each value's corresponding cell pairs
    v = [A,r,c];
end