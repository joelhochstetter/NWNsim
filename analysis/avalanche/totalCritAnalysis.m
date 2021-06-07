function totalCritAnalysis(baseFolder, subfolders, importMode, fitML, ncpu, whatMethods, Gthresh, relthresh, noisefloor, prename, binSize, conditions)
%{
Given a base folder and subfolders:
- Runs criticality analysis of each file from subfolder (either .mat,
    .tdms, .txt)
- Runs joint criticality analysis on each subfolder
- Runs joint criticality analysis on all data combined
- methods is a vector [1,1,1] means run all. [1,0,0] is just threshold
- conditions (struct): allows to split time-series pre-activation vs
    post activation, etc.
        e.g. to select post-activation times
        conditions.type   = 'crossing'
        conditions.after  = true
        conditions.thresh = 5e-6

%}

    %% Initialisation
    cd(baseFolder);
    
    if nargin < 10
        prename = '';
    end
    
    if nargin < 11
        binSize = -1;
    end
    
    if nargin < 12
        conditions = struct('type', 'none');
    end
    
    
    if ~iscell(subfolders)
        subfolders = {subfolders};
    end

    
    %% methods
    methods = {'threshold', 'ratioThreshold','hybrid', 'thresholdPeak', 'stationaryPt'};
    methods = methods(whatMethods);
    
    
    %% events
    eventDetect = struct();
    eventDetect.thresh     = Gthresh;
    eventDetect.relThresh  = relthresh;
    eventDetect.noiseFloor = noisefloor;       
    
    
    %% run
    for j = 1:numel(methods)
        for i = 1:numel(subfolders)
            eventDetect = struct();
            eventDetect.thresh     = Gthresh;
            eventDetect.relThresh  = relthresh;
            eventDetect.noiseFloor = noisefloor;                 
            
            eventDetect.method = methods{j};
            importFolder = strcat(baseFolder, '/', subfolders{i});
            saveFolder = strcat(prename, 'Avalanche_', num2str(methods{j}), '_Gt', num2str(eventDetect.thresh), '_rt', num2str(eventDetect.relThresh), '_nf', num2str(eventDetect.noiseFloor));
            multiCritAnalysis(importFolder, saveFolder, importMode, eventDetect, fitML, binSize, ncpu, conditions)
            cd(importFolder)
            saveFolder = strcat(prename, 'Avalanche_Joint_', num2str(methods{j}), '_Gt', num2str(eventDetect.thresh), '_rt', num2str(eventDetect.relThresh), '_nf', num2str(eventDetect.noiseFloor));        
            jointCritAnalysis(importFolder, saveFolder, importMode, eventDetect, fitML, binSize, conditions);
        end
        
        cd(baseFolder)
        saveFolder = strcat(prename, 'Avalanche_Joint_', num2str(methods{j}), '_Gt', num2str(eventDetect.thresh), '_rt', num2str(eventDetect.relThresh), '_nf', num2str(eventDetect.noiseFloor));
        jointCritAnalysis(subfolders, saveFolder, importMode, eventDetect, fitML, binSize, conditions);
    end




end