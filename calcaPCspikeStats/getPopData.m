function popData = getPopData(data)
    popData = cell(3,1);
    for i = 1:length(popData)
        for j = 1:length(data{i})
            popData{i} = vertcat(popData{i},transpose(data{i}{j}));
        end
    end
end