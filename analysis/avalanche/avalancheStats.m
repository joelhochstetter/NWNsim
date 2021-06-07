function [sizeAv, lifeAv, avTime, branchAv] = avalancheStats(events, t, joinperiod)
%{
    Input:
        events (Nx1 array) - number of events at given time bin

    Avalanche defined such that events happen in subsequent time bins and
    no events happen in preceding and after time-bin

    A avalanches are recorded

    If no timevector specified assume if it is 1:numel(events)
    
    joinperiod: for time series which join an ensemble of
    different simulations stores periodicity so ignores events calculated 
    between adjacent simulations


    Output:
         sizeAv (Ax1 array) - number of events in given avalanche
           lifeAv (Ax1 array) - length of avalanche (number of bins)
        avTime (Ax1 array) - time-stamp of the start of avalanche. Last
            blank frame before avalanche starts
     branchAv (Ax1 array) - (num events in time bin 2)/(num events in time bin 1)
       
%}


    if nargin < 2
        t = 1:numel(events);
    end

    if nargin < 3
        joinperiod = -1;
    end
    
    if joinperiod == -1
        joinperiod = numel(events) + 1;
    end    
    
    
    %handle joinperiod as a vector:
    if numel(joinperiod) == 1
        joinperiod = [0; [joinperiod:joinperiod:(numel(events) + 1)]'];
    end
    
    if joinperiod(1) > 1
        joinperiod = reshape(joinperiod, [numel(joinperiod), 1]);
        joinperiod = [0; joinperiod];
    end
    
    avEdg  = find(events == 0); %edges for avalanches
    
    A               = numel(avEdg) - 1;    
    sizeAv      = zeros(A,1);
    lifeAv        = zeros(A,1);
    branchAv = zeros(A,1);
    avTime     = zeros(A,1);
        
    joinID = discretize(avEdg - 1, joinperiod);
    
    for avId = 1:A
        % check they fall within the same join period
        % check that it's not an empty avalanche
        if joinID(avId) == joinID(avId + 1) && events(avEdg(avId) + 1) > 0
            sizeAv(avId)      = sum(events(avEdg(avId):avEdg(avId+1)));
            lifeAv(avId)        = avEdg(avId+1) - avEdg(avId) - 1;
            avTime(avId)     = t(avEdg(avId));
            branchAv(avId) = events(avEdg(avId)+2)/events(avEdg(avId)+1);                
        else
            sizeAv(avId)     = 0;
            lifeAv(avId)       = 0;    
            avTime(avId)    = 0;
            branchAv(avId) = 0;                
        end
    end   
    
    include     = sizeAv > 0;
    sizeAv      = sizeAv(include);
    lifeAv        = lifeAv(include);   
    avTime     = avTime(include);        
    branchAv = branchAv(include);    

end