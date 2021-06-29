function [G, V, t, I] = applyConditions(G, V, t, conditions)
%{
    Applies conditions as specified by the conditions struct to the
    conductance (G), voltage (V) and time (t) time-series
    
    e.g. 
    conditions = struct('type','crossing','after',true,'thresh'1e-6)
    conditions = struct('type','tInterval','lc', 1, 'uc', 5)
    
    Another option is allows decimation of data to a specified
        dt or to a certain factor
%}

    G = reshape(G, [numel(G),1]);
    t  = reshape(t, [numel(t),1]);
    V = reshape(V, [numel(V),1]);    
    
    if ~isfield(conditions, 'type')
        conditions.type = 'none';
    end
 
    if ~isfield(conditions, 'halfopen') %ends at last time-point
        conditions.halfopen = false;
    end    
    
    %another option is to decimate data
    %must supply either a 'new_dt' or destroyFactor
    if isfield(conditions, 'decimate') && (isfield(conditions, 'new_dt') || isfield(conditions, 'destroyFactor'))
        if isfield(conditions, 'new_dt')
            dt = (t(end) - t(1))/(numel(t) - 1);
            conditions.destroyFactor = round(conditions.new_dt/dt);
        end
        if conditions.destroyFactor > 1
            if conditions.decimate == 1 %decimate by skipping
                G = G(1:conditions.destroyFactor:numel(G));
                t  = t(1:conditions.destroyFactor:numel(t));            
            elseif conditions.decimate == 2 %decimate by averaging
                [G, t] = binData(G, conditions.destroyFactor);
                G = G/conditions.destroyFactor;
            elseif conditions.decimate == 3 %decimate with low-pass folder
                G = decimate(G, conditions.destroyFactor);      
                t  = t(1:conditions.destroyFactor:numel(t));            
            end
        end
    end
    
    I = 1:numel(G);
    %apply conditions for extracting an interval from simulations
    switch conditions.type
        case 'crossing'
            I = subTimeseries(G, conditions.after, conditions.thresh);
        case 'tInterval'
            I = extractInterval(t, conditions.lc, conditions.uc, false);
        case 'GInterval'
            I = extractInterval(G, conditions.lc, conditions.uc, false);  
        case 'Gchange'
            I = changeInterval(G, conditions.thresh);
        case 'dGGchange'
            I = eventInterval(G, inf, conditions.thresh, false);               
        case 'eventInterval'
             I = eventInterval(G, conditions.thresh, conditions.ratio, conditions.halfopen);
    end

    G = G(I);
    if numel(V) > 1
        V = V(I);
    end
    t = t(I);

end