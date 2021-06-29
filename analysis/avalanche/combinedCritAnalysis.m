function combinedCritAnalysis(netC, voltage, times, saveFolder, eventDetect, fitML, binSize, conditions)
%{
    Performs criticality analysis on
    
    Inputs:
        netC (cell array): each cell corresponds to network conductance
                time-series
        voltage (array): DC voltage for each time-series
        times(cell array): each cell corresponds to the time-vector for
            each time-series
        saveFolder (string): name of where to save.  e.g. 'avalancheAnalysis'
        eventDetect (how to do the event detection): struct:
            bareThreshold, dG./G threshold, peaks of dG/dt with appropriate
            smoothing, peaks with a 'refractory period'
            'method' = 'threshold', 'ratioThreshold', 'stationaryPt'
            'window' = running a running window mean
        conditions (struct): allows to split time-series pre-activation vs
            post activation, etc.
                e.g. to select post-activation times
                conditions.type   = 'crossing'
                conditions.after  = true
                conditions.thresh = 5e-6


    Written by Joel Hochstetter
%}         


    %% Set-up
    mkdir(saveFolder)
    
    if nargin < 8
        conditions = struct('type', 'none');
    end
    
    
    %% Process files and extract G, V, t
    Gjoin = [];
    tjoin = [0];
    Vjoin = [];
    fname = 'joint folders:';
    joinSpots = [0];
    
    saveNetC = true;    
    
    assert(numel(netC) == numel(times))
    
    for i = 1:numel(netC)
        assert(numel(netC{i}) == numel(times{i}));        
        [G, V, t]      = applyConditions(netC{i}, voltage(i), times{i}, conditions);
            
        if numel(t) == 0
            continue
        end

        t  = reshape(t, [1, numel(t)]);
        G = reshape(G,  [1, numel(G)]);
        V = reshape(V,  [1, numel(V)]);

        t = t + tjoin(end) - t(1);


        tjoin = [tjoin, t];
        Gjoin = [Gjoin, G];
        Vjoin = [Vjoin, V]; 
        joinSpots = [joinSpots, numel(t) + joinSpots(end)]; %stores index of joins (index of last element in time-series) for fixing IEI        
    end
    
    %cut off initial zero
    tjoin = tjoin(2:end);
    joinSpots = joinSpots(2:end);
    
    % detect events
    events =  findEvents(Gjoin, eventDetect);
    
    if numel(tjoin) == 0
       disp('No data entered') 
       return
    end
    
    dt = (tjoin(end) - tjoin(1))/(numel(tjoin) - 1);
    
    %perform criticality analysis
    critAnalysis(events, dt, Gjoin, tjoin, Vjoin, fname, strcat(saveFolder, '/'), fitML, binSize, joinSpots, saveNetC);

    
end