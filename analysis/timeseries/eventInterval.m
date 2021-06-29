function  [I, cuts] = eventInterval(x, thresh, ratio, halfopen)
%{
    Extracts the interval where (dx > thresh) || (dx/x > ratio) first time 
    to  (dx > ratio) || (dx/x > thresh) the last time. dx = gradient(x)


that occurs in a given time-series vector
    x based on a lower and upper cut-off 

    Inputs:
        x    (array): time-series vector
       thresh (double): threshold
       ratio     (double): ratio threshold on x
       halfopen (boolean): false: extract interval, true: end at the
            end of the time-series

    Outputs:
        I        (array): Interval containing the specified timepoints
        cuts (1x2 array): First and last index of time point 
%}

    I         = [];
    cuts = [nan, nan];

    assert(sum(isnan(x)) == 0);
    dx = gradient(x);
    pdx = find(abs(dx) > thresh | abs(dx./x) > ratio);
    if sum(pdx) > 0
        cuts(1) = pdx(1);
        if halfopen 
            cuts(2) = numel(x);
        else
            cuts(2) = pdx(end);    
        end
        I = cuts(1):cuts(2);
    end    

    
end