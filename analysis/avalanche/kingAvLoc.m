function kingLoc = kingAvLoc(x, N, c)
%{
    Determines location of king avalanche given bin-centre x and N (which can
    be either counts or probability.

    Method is find first increasing region. Then find the maxima over this
    region

    Inputs:
        x: bin centre of histogram
        N: Counts or probability assigned to each bin
        c: lower cut-off for which kingLoc can take. i.e. kingLoc > c 

    Outputs:
        kingLoc: bin centre corresponding to kind avalanche


    Written by Joel Hochstetter
%}

    if nargin < 3
        c = 0;
    end
    
    c = max([find(diff(N) > 0, 1), c]);
    N = N(x > c);
    x = x(x > c);
    
    if numel(N) > 0
        [~, idx] = max(N);
        kingLoc = x(idx);
    else
        kingLoc = nan;
    end
end