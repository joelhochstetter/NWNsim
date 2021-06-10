function [ieiData, ieiTime] = IEI(events, dt, joinperiod, t)
%{
    Calculates inter-event interval given  binarised
    time-series data and time-step
   ieiTime associated with  first event from IEI pair
   joinperiod: for time series which join an ensemble of
    different simulations stores periodicity so ignores events calculated 
    between adjacent simulations. Send in as an array if not periodic

    runMode = 1: use time-step
    runMode = 2: use time-vector to calculate IEI


    Written by Joel Hochstetter
%}

    if nargin < 4
        runMode = 1; %
    else
        runMode = 2;
    end
    
    if nargin < 3
        joinperiod = -1;        
    end
    
    if joinperiod == -1
        joinperiod = numel(events) + 1;
    end
        
    ieiData   = [];
    ieiTime  = []; %event times
    ieiIdx    = 1;
    prevEvent = find(events, 1);
    
    %joinperiod
    if numel(joinperiod) == 1
        joinperiod = joinperiod:joinperiod:(numel(events) + 1);
    end
    
    j = 1; %index for joinperiod
    
    if runMode == 1
        for i = (prevEvent + 1):numel(events)
            if events(i)
                if floor((i-1)/joinperiod(j)) == floor((prevEvent-1)/joinperiod(j))
                    ieiData(ieiIdx) = i - prevEvent;
                    ieiIdx = ieiIdx + 1;
                    ieiTime(ieiIdx) = prevEvent;
                    %increment j if above the join period
                    if ((i-1) > joinperiod(j)) && ((prevEvent-1) > joinperiod(j)) 
                        j = j + 1;
                    end
                end
                prevEvent = i;
            end
        end

        ieiData = ieiData*dt;
    elseif runMode == 2
        for i = (prevEvent + 1):numel(events)
            if events(i)
                if floor((i-1)/joinperiod(j)) == floor((prevEvent-1)/joinperiod(j))
                    ieiData(ieiIdx)  = i - prevEvent;
                    ieiTime(ieiIdx) = prevEvent; 
                    ieiIdx = ieiIdx + 1;
                    %increment j if above the join period
                    if ((i-1) > joinperiod(j)) && ((prevEvent-1) > joinperiod(j)) 
                        j = j + 1;
                    end                    
                end
                prevEvent = i;                
            end
        end     
        ieiTime = ieiTime(ieiData > 0);        
        ieiData  = ieiData(ieiData > 0);
        ieiTime = ieiTime(ieiTime > 0);
        assert(numel(ieiTime) == numel(ieiData));
    end

end