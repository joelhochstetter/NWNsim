function maxes = seqMaxCounts(seq)
%{
    For a cell array containing sequence of variables at
    each time point. Calculate the max of the variable 
%}

    maxes = zeros(size(seq));

    for i = 1:numel(seq)
        maxes(i) = max([seq{i},0]);
    end
    
    maxes(isnan(maxes)) = 0;

end