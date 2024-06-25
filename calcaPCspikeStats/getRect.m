%% function 'getRect' takes in matrix 'A' and the dimensions of the subgroups (nOb and nPc) within it. 
%  Captures subscripts of elements in nOb x nPc (upper right rectangle) as
%  'r' and 'c' and returns 3 column matrix 'v', with col1 as the value
%  corresponding to subscripts 'r' and 'c'.

function v = getRect(A,nOb,nPc)
    % 'v' contains upper right rectangular matrix values as 'A', and
    %   each value's corresponding cell pairs
    % 'r' contains row indices of each element in nOb x nPc
    r = reshape(repmat(1:nOb,nPc,1),nOb*nPc,1);
    % 'c' contains column indices of each element in nOb x nPc
    c = reshape(repmat((1:nPc)',nOb,1),nOb*nPc,1);
    % use 'r' and 'c' to get corresponding elements
    rect = A(1:nOb,nOb+1:end);
    idx = sub2ind(size(rect),r,c);
    rect = rect(idx);
    v = [rect,r,c];
end