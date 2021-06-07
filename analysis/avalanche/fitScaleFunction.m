function [acoeff, gamma] = fitScaleFunction(dur, re_tm, scale)
%{
    Input: 
        dur:   lifetime of avalanches
        re_tm: rescaled time
        scale: scaling function


    Ouput:
        fit the function scaleFunction for a,b,c,d,e, gamma for all of the
        datasets. Merge datasets for each T
%}


    %% Convert datapoints into three 1-D arrays
    %Get size
    sz = sum(dur);  
    T1 = zeros(sz,1);
    t1 = zeros(sz,1);
    s1 = zeros(sz,1);
    
    %Convert to 1D arrays
    k = 1;
    for i = 1:numel(dur)
        for j = 1:numel(scale{i})
            T1(k) = dur(i);
            t1(k) = re_tm{i}(j);
            s1(k) = scale{i}(j);
            k = k + 1;
        end
    end
    
    %% Calculates the sum of square of errors
    % a(6) = gamma
    err = @(a) sum((scaleFunction(t1, T1, a(1:5), a(6)) - s1).^2);
    a = fminsearch(err, [-1,0,0,0,0,1]);    
    acoeff = a(1:5);
    gamma  = a(6);
    
end