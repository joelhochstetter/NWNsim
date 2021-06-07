function [binned, times] = binData(data, tStep, roundUp)
%{
    Puts data into bins according with the tStep an integer corresponding
    to how many time-points per bin
    roundUp: is whether to round up number of bins. 
        default = false. Round down so exclude elements 
        in the last few bins

    Outputs:
            binned: sum of the data that goes in that bin
              times: times for the centres of each bin as an index
%}

    if nargin < 3
        roundUp = false;
    end

    len   = numel(data);
    nStep = floor(len/tStep);
    
    binned = zeros(nStep, 1);
    
    for i = 1:nStep
        binned(i) = sum(data((1 + (i - 1)*tStep):i*tStep));
    end

    times = [round(tStep/2):tStep:nStep*tStep]';
    
    % A bin of smaller size is added if number of bins does
    %not divide length and specify to
    if mod(len,tStep) > 0 && roundUp
        binned(nStep+1) = sum(data((1 + nStep*tStep):end));
        te = time(nStep) + tStep;
        if te > len
            times(nStep+1) = len;
        else
            times(nStep + 1) = te;
        end
    end
    
    assert(numel(times) == numel(binned));
    
end