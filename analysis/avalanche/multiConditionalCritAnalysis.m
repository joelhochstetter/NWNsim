function multiConditionalCritAnalysis(conditions, baseFolder, saveFolder, vals, subtype, binSize, fitML, joinperiod, fmt, dt)
%{
    For a set of simulations (note use jointCritAnalysis /
    multiCritAnalysis for experiments).

    multiConditionalCritAnalysis(struct('type','tInterval','lc', 2, 'uc', 8), '/import/silo2/joelh/Criticality/Avalanche/BigNetwork/RectElectChangeV/avNew', 'AvLikePhase', 1.05, 'Vstar', -1, false)

%}
    
    %% defaults   
    if nargin < 7
        fitML = false;
    end
    
    if nargin < 8
        joinperiod = 30*1000;
    end    
    
    if nargin < 9
        fmt = '%g';
    end

    if nargin < 10
        dt = 1e-3;
    end
    
    vals       = reshape(vals,       [1, numel(vals)]);
    binSize = reshape(binSize, [1, numel(binSize)]);    
    
    %%
    for v = vals
        for bs = binSize
            load(strcat(baseFolder, '/', subtype, num2str(v, fmt), '/events.mat'), 'events')
            load(strcat(baseFolder, '/', subtype, num2str(v, fmt), '/netC.mat'), 'netC')            
            saveFolder1 = strcat(baseFolder, '/', saveFolder, '/', subtype, num2str(v, fmt), '/bs', num2str(bs));            
            time = dt*(1:numel(events))';
            conditionalCritAnalysis(conditions, events, dt, netC, time, v, '', saveFolder1, fitML, bs, joinperiod);
        end
    end
end