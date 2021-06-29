function [bins, probs, edges] = LogBin(x, nbins, mx, R)
%{
    E.g:
        x = [1,1,1,1,1,5,10,12];
        LogBin(x, 2)
%}
    N = numel(x);
    x = reshape(x, [numel(x), 1]);
    %bins begin floor(c*R^j)
    c = 1;

    if nargin == 1
        R = 2;
        nbins = log(max(x))/log(c*R);        
    elseif nargin < 4
        x = reshape(x, [numel(x), 1]);
        %bins begin floor(c*R^j)
        if nargin == 2
            R = (max(x)^(1/nbins)/c);
        else
            R = (mx^(1/nbins)/c);
        end
        if R <= 1
            R = 2;
        end
    else
%         R = round(max(x)^(1/nbins)/c);
        
        nbins = log(max(x))/log(c*R);
    end
    
    
%     R = 2;
    edges = floor(c*R.^[0:(nbins + 1)]);
    edges = unique(edges);

%     bins = (edges(1:end - 1) + edges(2:end))/2;
%     bins1 = (edges(1:end - 1) + edges(2:end))/2;
    
    bins = floor(sqrt(edges(1:end - 1).*edges(2:end)));
%     bins  = floor(c*R.^[0.5:(nbins)]);
    bins    = unique(bins);
    bsize = edges(2:end) - edges(1:end -1);
%     bins(bsize == 1) = bins1(bsize == 1);
    probs = sum((x < edges(2:end)) & (x >= edges(1:end-1)))./bsize/N;


end
