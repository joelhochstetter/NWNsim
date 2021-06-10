function [binned,  times,  binJoinPeriod] = binWithJoins(data, tStep, joinperiod)
%{
    Puts data into bins according with the tStep an integer corresponding
    to how many time-points per bin
    roundUp: is whether to round up number of bins. 
        default = false. Round down so exclude elements 
        in the last few bins

    Outputs:
              binned: sum of the data that goes in that bin
                times: the central time-point for each bin
   binJoinPeriod: array of the join

    Written by Joel Hochstetter

%}
    

    if numel(joinperiod) <= 0
        binJoinPeriod = -1;
        binned = binData(data, tStep);
        return;
    elseif numel(joinperiod) == 1
        if joinperiod <= 0
            joinperiod = numel(data);
        end
        joinperiod = [0, joinperiod:joinperiod:numel(data)];
    elseif joinperiod(1) > 1
        joinperiod = [0, joinperiod];
    end
    
    data = reshape(data, [numel(data),1]);
    
    numPeriods = numel(joinperiod) - 1;
    
    binned = zeros(numPeriods*tStep, 1);
    times   = zeros(numPeriods*tStep, 1); 
    sidx = 1; %start index
    
    binJoinPeriod = zeros(numPeriods + 1,1);
    
    for j = 1:numPeriods
        N = floor((joinperiod(j+1) - joinperiod(j))/tStep) - 1;
        binJoinPeriod(j + 1) = sidx + N;
        [binned(sidx:sidx+N), times(sidx:sidx+N)] = binData(data(1+joinperiod(j):joinperiod(j+1)), tStep, false);
        times(sidx:sidx+N) = times(sidx:sidx+N) + joinperiod(j);
        sidx = sidx + N + 1;
    end
    
    assert(numel(times) == numel(binned))

end