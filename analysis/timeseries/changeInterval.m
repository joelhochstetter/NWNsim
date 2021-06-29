function [I, cuts] =  changeInterval(x,thresh)
%{
    Extracts the interval where dx > thresh first time 
    to dx > thresh the last time. dx = gradient(x)


that occurs in a given time-series vector
    x based on a lower and upper cut-off 

    Inputs:
        x    (array): time-series vector
       thresh (double): threshold

    Outputs:
        I        (array): Interval containing the specified timepoints
        cuts (1x2 array): First and last index of time point 
%}
    
    I         = [];
    cuts = [nan, nan];

    assert(sum(isnan(x)) == 0);
    dx = gradient(x);
    pdx = find(abs(dx) > thresh);
    if sum(pdx) > 0
        cuts(1) = pdx(1);
        cuts(2) = pdx(end);    
        I = cuts(1):cuts(2);
    end    
end