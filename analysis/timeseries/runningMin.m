function mx = runningMin(x, dt)
%{
    Given a time-series calculates a running minimum, bins into central
    time-point of time series (if dt is odd) and one before if dt is even
    %}

    %% Initialise outputs
    mx = zeros(size(x));    
    N  = numel(x);
    dt1 = floor((dt - 1)/2); %lower range
    dt2 = dt - dt1 - 1;      %upper range
    
    %% Initial datapoints
    for i = 1:dt1
        mx(i) = min(x(1:i));
    end   
    
    %% Central datapoints
    for i = (dt1 + 1):(N-dt2)
        mx(i) = min(x(i - dt1:i + dt2));
    end
    
    %% End datapoints
    for i = (N - dt2 + 1):N
        mx(i) = min(x(i:end));
    end     


end