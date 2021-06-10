function jointCritAnalysis(importFolder, saveFolder, importMode, eventDetect, fitML, binSize, conditions)
%{
    Looks through data and performs a criticality analysis on all files in
    the folder. This works for experimental data and simulated data
    
    Inputs:
        importFolder (string) or cell: folder to extract data from. If
            enter a cell then loops over all folders
        saveFolder (string): name of where to save.  e.g. 'avalancheAnalysis'
        importMode (int): 0 = simulation file, 1 = TDMS file, 2 = text file   
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
    
    if ~iscell(importFolder)
        importFolder = {importFolder};
    end
    
    if nargin < 7
        conditions = struct('type', 'none');
    end
    
    %% Process files and extract G, V, t
    Gjoin = [];
    tjoin = [0];
    Vjoin = [];
    fname = 'joint folders:';
    joinSpots = [0];
    
    saveNetC = true;    
    
    numFolders = numel(importFolder);
    for j = 1:numFolders
        numFiles = howManyFiles(importMode, importFolder{j});
        fname = strcat(fname, importFolder{j}, ',');
        for i = 1:numFiles
            %import file
            [G, V, t, ~] = importByType(importMode, importFolder{j}, i);
            [G, V, t]    = applyConditions(G, V, t, conditions);
            
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
    end
    
    %cut off initial zero
    tjoin = tjoin(2:end);
    joinSpots = joinSpots(2:end);
    
    % detect events
    events =  findEvents(Gjoin, eventDetect);
    if numel(tjoin) == 0
       disp('None') 
       return
    end
    dt = (tjoin(end) - tjoin(1))/(numel(tjoin) - 1);
    
    %perform criticality analysis
    critAnalysis(events, dt, Gjoin, tjoin, Vjoin, fname, strcat(saveFolder, '/'), fitML, binSize, joinSpots, saveNetC);

    
end