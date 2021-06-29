function var = runningVariability(x, dt, alpha)
%{
    Calculates a measure of variability across different parts
    of the time-series
    
    alpha = 1: Coefficient of variation (sd^1/mean)
    alpha = 2: Fano factor (sd^2/mean) 
    
    Inputs:
        x:  time-series data, normally conductance
        dt: number of time-points which we group together
     alpha: the exponent used to determine the variablity quantity
     calculated
    
    Output:
        res: the running exponents for each time-point. Start and end
        time-points use first dt time-steps.
%}
    %% Initialise outputs
    var = zeros(size(x));  %variability coefficient  
    N  = numel(x);
    dt1 = floor((dt - 1)/2); %lower range
    dt2 = dt - dt1 - 1;      %upper range
    
    %% Initial datapoints
    var(1) = std(x(1:dt)).^alpha/mean(x(1:dt));
    for i = 2:dt1
        var(i) = var(1);
    end   
    
    %% Central datapoints
    for i = (dt1 + 1):(N-dt2)
        var(i) = std(x(i - dt1:i + dt2)).^alpha/mean(x(i - dt1:i + dt2));
    end
    
    %% End datapoints
    var(N) = std(x(end - dt + 1:end)).^alpha/mean(x(end - dt + 1:end));

    for i = (N - dt2 + 1):(N - 1)
        var(i) = var(N);
    end     


end