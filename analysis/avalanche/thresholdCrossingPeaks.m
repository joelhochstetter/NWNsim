function events = thresholdCrossingPeaks(x, thresh)
    events = zeros(size(x));
    [~, above, below] = thresholdCrossing(x, thresh);
    %assumes time-series does not start above threshold else cuts
    %to this segement
    if numel(above) > 2 && numel(below) > 2
        if above(1) > below(1)
            above = above(2:end);
            below = below(1:end-1);
        end
        if numel(above) > numel(below)
           above = above(1:end-1); 
        elseif numel(above) < numel(below)
           below = below(1:end-1); 
        end
        assert(above(1) < below(1));
        assert(numel(above) == numel(below))
        for i = 1:numel(above)
            [~, I] = max(x(above(i):below(i)));
            events(I + above(i) - 1) = true;
        end
    end    

end