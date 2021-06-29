function DC_Vsweep(saveFolder, Amps, T, connFile, initStateFile , initStateFolder, initCon, rescalePLength, useRect, rectFraction, contactDistance, saveFilState, saveEvents, EventThreshold, nameComment)
%{
    Runs DC sweep over different voltages

    Compulsory inputs:
        saveFolder (string): Folder to save simulations
        Amps         (Nx1 dbl): DC amplitude for voltage sweep in V
        T               (double): Length of time to simulate for


    Optional inputs:
        Network connectivity:
            connFile   (string): name of file containing connectivity, defaults
                to 100 nanowire, 261 junction 30umx30um

        Initial junction states:            
            initStateFile (string): Name of file to which to use initial
                states. If want initial states all zero, then initStateFile = 0
            initStateFolder (string): If providing file of initial states.
                Defaults to ".".  Only required if initStateFile provided
            initCon (double): If not provided or set to -1: then we take the initial states
                from the end of the simulation given by initStateFile
                If provided: then chooses the junction states from the
                time-point where conductance is closest to initCon

        Electrodes used:
            rescalePLength (double): re-scale stimulus voltage by
                source-drain path length. So V = V*sdPathLength. default: 0
           useRect (boolean): false => use point electrodes, true => use
                rectangular electrodes at opposite ends of NWN. Defaults to false
            rectFraction (double): vertical fraction of network that
                rectangular electrode covers. defaults to 0.035
            contactDistance (integer): for point electrode specifies the
                distance between electrodes as a source-drain (topological)
                path length. by default (or if contactDistance < 0) 
                takes physically furthest path length.


        What to save:
            saveFilState (double): Specify whether or not to save junction
                filament states at each time-step. False saves at the final
                time-point only
            saveEvents (boolean): true: saves events of |dG/dt| exceeding
                threshold before falling below threshold
            EventThreshold (double): positive number specifies the r by
                which |dG/dt| > r before returning below r defines an event
            nameComment (string): string to append to end of the name of '.mat' file

    Outputs:
        Runs simulations for each of specified voltages and outputs
            simulation files into saveFolder


    Written by Joel Hochstetter

%}

    %% Set-up defaults for function arguments.
    %must provide first 3 args (saveFolder, Amps, T)
    
    
    %default nameComment
    if nargin < 15
        nameComment = '';
    end
    
    %initial connectivity
    if nargin < 4 ||  (isnumeric(connFile) && connFile <= 0)
        connFile = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
    end
    
    %specify closest initial conductance value
    %if -1: defaults to end of simulation if initStateFile provided
    if nargin < 7
        initCon = -1;
    end
    
    if nargin < 5 || (isnumeric(initStateFile) && initStateFile <= 0)
        initLamda = 0;
    else
        if nargin == 5
            initStateFolder = saveFolder;
        end
        sim = multiImport(struct('SimOpt', struct('saveFolder', initStateFolder), 'importByName', initStateFile, 'importStateOnly', true));
        if isfield(sim{1}, 'swLam')
            if initCon < 0
                initLamda = sim{1}.swLam(end,:)';
                nameComment = strcat(nameComment, '_initCon', num2str(sim{1}.netC(end)), 'S');
            else
                [~, tidx] = min(abs(sim{1}.netC - initCon));
                initLamda = sim{1}.swLam(tidx,:)';
                nameComment = strcat(nameComment, '_initCon', num2str(sim{1}.netC(tidx)), 'S');
            end
        elseif isfield(sim{1}, 'finalStates')
            initLamda                  = sim{1}.finalStates';
        else
            disp('FAILED');
            return;
        end
    end
    
    if nargin < 8
        rescalePLength = false;
    end

    if nargin < 9
        useRect = 0;
    end

    if nargin < 10
        rectFraction = 0.035;
    end
    
    if nargin < 11 || contactDistance < 0
        contactMode =  'farthest';
        contactDistance = -1;
    else
        contactMode = 'fixedTopoDistance';
    end
    
    if nargin < 12
        saveFilState = true;
    else
        saveFilState = boolean(saveFilState);
    end

    if nargin <13
        saveEvents = true;
    else 
        saveEvents = boolean(saveEvents);
    end

    if nargin < 14
        EventThreshold = 1e-3;
    end


    %% Set-up parameters for simulation
    params = struct();

    % Set Simulation Options
    params.SimOpt.useWorkspace    = false; %only save simulations do need keep values in workspace
    params.SimOpt.saveSim         = true;
    params.SimOpt.useParallel     = true; %to use a parallel simulations specify this as true
    %if you want to run only one simulation at once then you can specify idx
        %and uncomment next line
%     params.SimOpt.runIndex = idx; 
    params.SimOpt.hdfSave         = ~saveEvents & saveFilState; %If not saving events and saving filament state use hdf5 file format
    params.SimOpt.saveSwitches = false;

    params.SimOpt.saveFilStateOnly = saveFilState;
    params.SimOpt.saveEventsOnly  = saveEvents;
    params.SimOpt.EventThreshold  = EventThreshold;

    params.SimOpt.stopIfDupName = true; %this parameter only runs simulation if the savename is not used.
    params.SimOpt.saveFolder      = saveFolder;
    mkdir(params.SimOpt.saveFolder);

    params.SimOpt.T                = T;
    params.SimOpt.dt               = 1e-3;
    params.SimOpt.ContactMode = contactMode;
    params.SimOpt.ContactGraphDist = contactDistance;

    %Set Stimulus
    params.Stim.BiasType     = 'DC'; 
    params.Stim.Amplitude    = Amps; 

    %Set connect file
    params.Conn.filename = connFile;

    % use rectangular electrodes instead of point electrodes
    if useRect
        %Connectivity and contacts
        Connectivity          = getConnectivity(params.Conn);
        params.SimOpt.RectElectrodes = true;
        params.SimOpt.NewEdgeRS      = false; %edges connecting electrodes to nanowires are not memristive
        params.SimOpt.RectFractions  = rectFraction;
        params.SimOpt.XRectFraction  = 1.0;    
        [~, ~, SDpath] = addRectElectrode(Connectivity, params.SimOpt.RectFractions, params.SimOpt.XRectFraction);
        %Checks that a path containing
        if SDpath == Inf 
            disp('No such SD path')
            return
        end
    end    

    %if re-scaling path lengths extract source-drain path lengths (SDpath) then
        %multipy Voltages by SDpath
    if rescalePLength
        if ~useRect
            [Connectivity] = getConnectivity(params.Conn);
            SimulationOptions = selectContacts(Connectivity, params.SimOpt);
            Contacts = SimulationOptions.ContactNodes;
            SDpath = distances(graph(Connectivity.weights), Contacts(1), Contacts(2));
        end
        params.Stim.Amplitude = SDpath*params.Stim.Amplitude;
    end

    params.SimOpt.nameComment = nameComment;

    %Set Components
    params.Comp.ComponentType  = 'tunnelSwitchL';
    params.Comp.onConductance   = 7.77e-5;
    params.Comp.offConductance   = 7.77e-8;
    params.Comp.setVoltage     = 1e-2;
    params.Comp.resetVoltage  = 5e-3;
    params.Comp.criticalFlux   =  0.01;
    params.Comp.maxFlux        = 0.015;
    params.Comp.penalty        =    1;
    params.Comp.boost          =   10;
    params.Comp.filamentState = initLamda;


    %% Run simulations
    multiRun(params);


end


