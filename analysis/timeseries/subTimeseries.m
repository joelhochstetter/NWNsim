function [I, cuts] = subTimeseries(X, after, thresh)
%{
    extracts an interval of time-series based on the value of first time to
    cross a given threshold.

    Inputs:
        after  (binary): true: extracts interval after cut-off first crossed 
                        false: extracts interval before cut-off first
                            crossed
        thresh (double): threshold that time-series cross

    Outputs:
        I        (array): Interval containing the specified timepoints
        cuts (1x2 array): First and last index of time point 
%}
    I    = [];
    cuts = [];
    x0 = find(X >= thresh, 1);
    if numel(x0) > 0
        if after
           I    = x0:numel(X);
           cuts = [x0, numel(X)];
        else
           x0 = x0 - 1;
           if x0 > 0
               I    = 1:x0; 
               cuts = [1, x0];       
           end
        end
    end
end