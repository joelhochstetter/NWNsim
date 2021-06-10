function branch = branchingRatio(events, binSize)
%{
    Calculates the branching ratio given binSize
        
    As defined in Beggs Plenz J. Neuroscience (2003) 
        and de Carvalvo Prado PRL (2000)
    
    Written by Joel Hochstetter

%}

    binned = binData(events, binSize);
    [~, size_t, ~, numAv] = avalancheShape(binned);
    
    branches = zeros(size(size_t)); %size_t(1) is always 0
    for i = 1:numel(size_t)
        branches(i) = size_t{i}(3)/size_t{i}(2);
    end
    branch = sum(branches.*numAv)/sum(numAv);
    
end