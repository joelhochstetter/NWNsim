function conditionalCritAnalysis(conditions, events, dt, G, time, V, filename, saveFolder, fitML, binSize, joinperiod)
%{
    Run avalanche analysis limited by certain conditions such as time,
    conductance, etc.

%}

    %% Modify joinperiod so is in right format to process
    if joinperiod == -1
        joinperiod = numel(events) + 1;
    end
    
    if numel(joinperiod) == 1
        joinperiod = [0; [joinperiod:joinperiod:(numel(events) + 1)]'];
    end
    
    if joinperiod(1) > 1
        joinperiod = reshape(joinperiod, [numel(joinperiod), 1]);
        joinperiod = [0; joinperiod];
    end    

    %% Applying the conditions
    G1          = [];
    events1 = [];
    t1           = [];
    numPeriods = numel(joinperiod) - 1;
    conJoinPeriod = zeros(numPeriods + 1, 1);
    
    for i = 1:numPeriods
        [cG, ~, ct, I] = applyConditions(G(1+joinperiod(i):joinperiod(i+1)), 0, time(1+joinperiod(i):joinperiod(i+1)) - joinperiod(i)*dt, conditions);
        G1 = [G1; cG];
        t1   = [t1; ct];
        iEvents  = events(1+joinperiod(i):joinperiod(i+1));
        events1 = [events1; iEvents(I)];
        conJoinPeriod(i+1) = conJoinPeriod(i) + numel(cG);
    end
    
    
    %%
    critAnalysis(events1, dt, G1, t1, V, filename, saveFolder, fitML, binSize, conJoinPeriod, false);
    
    
end