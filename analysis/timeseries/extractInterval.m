function [I, cuts] =  extractInterval(x, lc, uc, long)
%{
    Extracts the first interval that occurs in a given time-series vector
    x based on a lower and upper cut-off 

    Inputs:
        x    (array): time-series vector
       lc (double): lower cut-off of variable x
      uc (double): upper cut-off of variable x
  long (double): extracts longest interval

    Outputs:
        I        (array): Interval containing the specified timepoints
        cuts (1x2 array): First and last index of time point 
%}
    
    I         = [];
    cuts = [nan, nan];
    if sum(isnan(x)) > 0
        assert(sum(isnan(x)) == 0);
    end
    onInt    = x >= lc & x <= uc;
    if sum(onInt) > 0 %checks interval is non-empty
        cuts(1) = find(x >= lc & x <= uc, 1);
        onInt = onInt(cuts(1):end);
        if all(onInt) %check whether interaval extends to end of time-series
            cuts(2) = numel(x);
        else
            if long  %extract longest interval
                f =  find(~onInt);
                cuts(2) = f(end) + cuts(1) - 2; %find last element off interval after onInt
            else %extracts first interval
                cuts(2) = find(~onInt, 1) + cuts(1) - 2; %find 1st element off interval after onInt                
            end
        end
        I = cuts(1):cuts(2);
    end
    
end