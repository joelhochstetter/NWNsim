function [dur, size_t, time_t, numAv, sample, lifeAv] = avalancheShape(events, sampleFraction, seed)
%{
    Input:
        events (Nx1 array) - number of events at given time bin

    Avalanche defined such that events happen in subsequent time bins and
    no events happen in preceding and after time-bin

    A avalanches are recorded

    Can bootstrap by providing a sample fraction and seed


    Output:
        dur (Nx1 array) - each unique avalanche duration
        size_t (Nx1 cell) - average avalanche size as a function of time
                              for each unique duration. Each element is a
                              time-series vector
        time_t (Nx1 cell) - stores the time vectors for the duration of
                              each event in the analysis

        numAv (Nx1 array) - stores number of avalanches of this length.

    Written by Joel Hochstetter       
%}

    if nargin < 2
        sampleFraction = 1.0;
    end
    
    if nargin < 3
        seed = -1;
    end

    minNum = 0; %minimum number of avalanches

    avEdg  = find(events == 0); %edges for avalanches   
    A = numel(avEdg) - 1;
    lifeAv = zeros(A,1);
    
    for avId = 1:A
        lifeAv(avId) = avEdg(avId+1) - avEdg(avId) - 1;
    end   
   
    
    
    %% boostrap
    nzAv = find(lifeAv > 0); %indices of non-zero avalanches
    Anz = numel(nzAv); %number of actual avalanches
    
    if seed >= 0
        %number of avalanches sampled 
        AS = round(sampleFraction*Anz);

        %Generate samples
        rng(seed);    
        sample = randi(Anz, AS, 1);
        
        sampleFull = nzAv(sample); %includes zero avalanches
        lifeAv = lifeAv(sampleFull); 
        
    else
        sampleFull = 1:A; %sample all avalanches
        sample = 1:Anz;
    end
    
    
    %% generate durations
    dur = unique(lifeAv(lifeAv > 0), 'sorted');
    N   = numel(dur); 
    
    size_t = cell(N,1);
    time_t = cell(N,1);
    numAv  = zeros(N,1);
    
    %make time vectors
    for i = 1:N
        time_t{i} = 0:(dur(i) + 1); %duration + timesteps before and after
    end
    
    
    %%
    %make average event at each time point
    for i = 1:N %loop over event durations
        size_t{i} = zeros(dur(i) + 2,1);
        avIDs = find(lifeAv == dur(i)); %get the avalanche IDs
        avIDs = sampleFull(avIDs); %convert back to original ID      
        nRelv = numel(avIDs); %Number of relevant avalanches
        if nRelv > 0
            for j = 1:(dur(i) + 2) %loop over size of avalanche            
                for k = 1:nRelv %loop over the relevant avalanches
                    if avEdg(avIDs(k)) + j - 1 > numel(events)
                        break;
                    end
                    size_t{i}(j) = size_t{i}(j) + events(avEdg(avIDs(k)) + j - 1);
                end
            end
        end
        size_t{i} = size_t{i}/nRelv; %convert sum to average
        numAv(i)  = nRelv;
    end
    
    %%
    lifeAv = lifeAv(lifeAv > 0);
    
   


end