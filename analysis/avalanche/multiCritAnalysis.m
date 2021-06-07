function multiCritAnalysis(importFolder, saveFolder, importMode, eventDetect, fitML, binSize, ncpu, conditions)
%{
    Looks through data and performs a criticality analysis on all files in
    the folder. This works for experimental data and simulated data
    
    Inputs:
        importFolder (string): folder to extract data from
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

%}         

    %% Defaults
    if nargin < 7
       ncpu = 1; 
    end
    
    if nargin < 8
        conditions = struct('type', 'none');
    end
    
    saveNetC = true;
    
    %start parallel pool if not already running
    if isempty(gcp('nocreate'))
        parpool(ncpu)
    end

    %% Set-up
    cd(importFolder)
    mkdir(fullfile(saveFolder))
    
    
    %% Process files and extract G, V, t
    numFiles = howManyFiles(importMode, importFolder);
%     for i = 1:numFiles    
    for i = 1:numFiles
        %import file
        [G, V, t, fname] = importByType(importMode, importFolder, i);
        dt = (t(end) - t(1))/(numel(t) - 1);
        [G, V, t] = applyConditions(G, V, t, conditions);
        
        if numel(t) == 0
            continue
        end

        % detect events
        events =  findEvents(G, eventDetect);
        
        %perform criticality analysis
        critAnalysis(events, dt, G, t, V, fname, strcat(saveFolder, '/', fname, '/'), fitML, binSize, -1, saveNetC);
    end
       
   
end