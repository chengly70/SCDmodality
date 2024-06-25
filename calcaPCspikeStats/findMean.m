%% function 'findMean' takes in varargin and returns 3x1 cell array corresponding to each drug state Bic, Mus, ND.
%  Each drug state cell contains another cell array to store means for each
%  recording. Each recording then has a cell array to store the means over
%  the 10 trials.
function c = findMean(varargin)
    c = cell(length(varargin),1);
    c{length(varargin),1} = [];
    for i = 1:length(varargin)
        v = cell(length(varargin{i}),1);
        v{length(varargin{i}),1} = [];
        for j = 1:length(varargin{i})
            v{j,:} = mean(varargin{i}{j,:});
        end
        c{i,:} = v;
    end
end