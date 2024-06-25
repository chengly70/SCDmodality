function v = findNumCells(A)
    v = zeros([length(A),1]);
    for i = 1:length(A)
       v(i,1) = width(A{i,1});
    end
end