function r = aveMaxMin(x, f, dt)
%{
    Calculates average ratio between maximum and 
        minimum conductance for a periodic driving signal

    Inputs:
        x (Nx1 double): time-series (e.g. conductance)
        f (double): frequency of the periodic driving signal
      dt (double): Simulation time-step in seconds


    Outputs:    
        r (double): average ratio between maximum and min of variable x


    Caveats: 1 period should correspond to an integer number of time-steps
        to avoid rounding errors


    Written by Joel Hochstetter
%}

    

    N = numel(x); %number of time-steps
    T = dt*N; %length of time of simulation
    p = 1/f;  %period of signal
    tstepP = p/dt; %number of time-steps in period
    
    if rem(p , dt)  >= 1e-8
        disp('WARNING: non-integer number of time-steps in period')
    end
        
    numP = floor(T/p);
    rP = zeros(numP,1); %max(x)/min(x) for each period
    
    for i = 1:numP
        maxX = max(x(round(tstepP*(i - 1) + 1):round(tstepP*i)));
        minX= min(x(round(tstepP*(i - 1) + 1):round(tstepP*i)));
        rP(i) = maxX/minX;
    end

    r = mean(rP);

end