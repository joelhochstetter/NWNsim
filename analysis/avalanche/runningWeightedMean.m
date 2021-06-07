function [mx, fx] = runningWeightedMean(x, w, dt)
%{
    Given a time-series calculates a running mean, bins into central
    time-point of time series (if dt is odd) and one before if dt is even
    
    w is the weight. We renormalise mean by wi/sum(w). Where we sum
    over the time-bin. Must be same size as x

    also returns a fluctuation about the mean at each point in time-bin.
%}

    %% Initialise outputs
    mx = zeros(size(x));    
    N  = numel(x);
    dt1 = floor((dt - 1)/2); %lower range
    dt2 = dt - dt1 - 1;      %upper range
    
    %% Initial datapoints
    for i = 1:dt1
        mx(i) = sum(x(1:i).*w(1:i))/sum(w(1:i));
    end   
    
    %% Central datapoints
    for i = (dt1 + 1):(N-dt2)
        mx(i) = sum(x(i - dt1:i + dt2).*w(i - dt1:i + dt2))/sum(w(i - dt1:i + dt2));
    end
    
    %% End datapoints
    for i = (N - dt2 + 1):N
        mx(i) = sum(x(i:end).*w(i:end))/sum(w(i:end));
    end 
    
    mx(isnan(mx)) = 0;
    
    %%
    fx = x - mx;

end