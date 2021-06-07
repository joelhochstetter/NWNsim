function [bins, probs, edges] = LinearBin(x, nbins)
%{
    E.g:
        x = [1,1,1,1,1, 5, 10, 12];
        LinearBin(x, 2)
%}

    N = numel(x);
    x = reshape(x, [numel(x), 1]);

    R = floor(max(x)/(nbins));
    
    edges = 1:R:(max(x) + R);
    bins  = (edges(1:end-1) + edges(2:end))/2;
    probs = sum((x < edges(2:end)) & (x >= edges(1:end-1)))/R/N;


end