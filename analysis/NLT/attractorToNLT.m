function NLTres = attractorToNLT(attractorFolder, saveFolder, nds)
%{
    From a folder of attractors performs the non-linear transform task of
    triangular wave to sine, square, pi/2 phase shifted triangular and double
    frequency triangular wave
    
    Inputs:
        attractorFolder (string): name of folder of attractors
        saveFolder (string): Folder to save NLT task results
        nds (vectors of integers): nodes to train for NLT task
        
    Outputs:
        Saves with MSE, RNMSE, target, actual resultant vector for each target signal
        
        Saves result to struct called 'NLTres'
        
        
    Written by Joel Hochstetter
%}


    if nargin < 3
        nds = 1:100;
    end

    if nargin < 4
        saveName = 'NLTresults.mat';
    end

    addpath(genpath(attractorFolder));
    files = dir(strcat(attractorFolder, '/*.mat'));
    
    dt = 1e-3;
    T  = 80;
    
    numTSteps = round(T/dt);
    
    Freqs = zeros(numel(files), 1);
    Amps = zeros(numel(files), 1);
    
    sinMSE      = zeros(numel(files), 1);
    sinRNMSE = zeros(numel(files), 1);
    sinResult  = zeros(numel(files), numTSteps); 
    
    dblMSE      = zeros(numel(files), 1);
    dblRNMSE = zeros(numel(files), 1);
    dblResult  = zeros(numel(files), numTSteps); 

    phsMSE      = zeros(numel(files), 1);
    phsRNMSE = zeros(numel(files), 1);
    phsResult  = zeros(numel(files), numTSteps);     

    squMSE      = zeros(numel(files), 1);
    squRNMSE = zeros(numel(files), 1);
    squResult  = zeros(numel(files), numTSteps);         
    
    parfor j = 1:numel(files)
        sim = multiImport(struct('SimOpt', struct('saveFolder', attractorFolder), 'importByName', files(j).name));
        
        %% set-up parameters
        params = struct();
        
        % Set Simulation Options
        params.SimOpt.saveSim         = true;
        params.SimOpt.onlyGraphics    = true; %does not plot anything
        params.SimOpt.useParallel     = false;
        params.SimOpt.hdfSave         = false;
        params.SimOpt.stopIfDupName = true; %this parameter only runs simulation if the savename is not used.

        params.SimOpt.dt               = dt;
        params.SimOpt.T                = T;        
        params.SimOpt.nameComment      = '';
        params.SimOpt.saveFolder       = strcat(attractorFolder, '/', saveFolder);
        mkdir(params.SimOpt.saveFolder)
        
        params.Comp.ComponentType = 'tunnelSwitchL';
        params.Stim.BiasType = sim{1}.Stim.BiasType;
        params.Stim.Amplitude = sim{1}.Stim.Amplitude;    
        params.Stim.Frequency = sim{1}.Stim.Frequency; 
        Freqs(j) = params.Stim.Frequency ;
        Amps(j) = params.Stim.Amplitude ;
        
        params.Comp.setVoltage = sim{1}.Comp.setV;
        params.Comp.resetVoltage = sim{1}.Comp.resetV;
        params.Comp.criticalFlux = sim{1}.Comp.critFlux;
        params.Comp.maxFlux = sim{1}.Comp.maxFlux;
        params.Comp.boost = sim{1}.Comp.boost;
        params.Comp.penalty = sim{1}.Comp.pen;
        
        params.SimOpt.ContactMode      = 'preSet';
        params.SimOpt.ContactNodes  = sim{1}.ContactNodes;
        
        params.Conn.filename = sim{1}.ConnectFile;     
        if isfield(sim{1}, 'swLam')
            initLamda                  = sim{1}.swLam(end,:)';
        elseif isfield(sim{1}, 'finalStates')
            initLamda                  = sim{1}.finalStates';
        else
            disp('FAILED');
            continue;
        end
        
        params.Comp.filamentState = initLamda;

        
        %% Run simulation
        multiRun(params);
        sim = multiImport(params);
        sim = sim{1};
        contacts = sim.ContactNodes;
        timeVec = sim.Stim.TimeAxis;
        connectivity = struct('filename', sim.ConnectFile);
        connectivity = getConnectivity(connectivity);
        swV = sim.swV; %extract junction voltages
        nwV = zeros(size(swV,1), connectivity.NumberOfNodes);
        
        if isfield(sim, 'nwV')
            nwV = sim.nwV;
        else
            for i = 1:size(swV, 1)
                nwV(i, :) = getAbsoluteVoltage(swV(i,:), connectivity, contacts);
            end           
        end
        
        
        %% Set-up target
        Amp = params.Stim.Amplitude;
        Frq   = params.Stim.Frequency;
        SigType =  params.Stim.BiasType;

        %targets: square, sawtooth, 2f-sine, cosine
        targetPhase = 0;
   
        %% 2-frequency
        targetStim = struct('Amplitude', Amp, 'Frequency', 2*Frq, 'Phase', targetPhase, 'BiasType', SigType);
        targetSOpt = struct('T', sim.T, 'dt', sim.dt);
        target = getStimulus(targetStim, targetSOpt);

        %%%
        [weights, mse, rnmse, y] = NLT(target.Signal, nwV, nds);
  
        dblMSE(j)      = mse;
        dblRNMSE(j) = rnmse;
        dblResult(j,:)  = y;     

        
        %% shifted out of phase
        targetStim = struct('Amplitude', Amp, 'Frequency', Frq, 'Phase', pi/2, 'BiasType', SigType);
        targetSOpt = struct('T', sim.T, 'dt', sim.dt);
        target = getStimulus(targetStim, targetSOpt);

        %%%
        [weights, mse, rnmse, y] = NLT(target.Signal, nwV, nds);
  
        phsMSE(j)      = mse;
        phsRNMSE(j) = rnmse;
        phsResult(j,:)  = y;             
        
        
        %% sine wave
        targetStim = struct('Amplitude', Amp, 'Frequency', Frq, 'Phase', targetPhase, 'BiasType', 'AC');
        targetSOpt = struct('T', sim.T, 'dt', sim.dt);
        target = getStimulus(targetStim, targetSOpt);

        %%%
        [weights, mse, rnmse, y] = NLT(target.Signal, nwV, nds);
 
        sinMSE(j)      = mse;
        sinRNMSE(j) = rnmse;
        sinResult(j,:)  = y;             

        %% square wave
        targetStim = struct('AmplitudeOn', Amp, 'AmplitudeOff', -Amp, 'OffTime', 1/(2*Frq), 'Phase', targetPhase, 'BiasType', 'Square', 'Duty', 50.0);
        targetSOpt = struct('T', sim.T, 'dt', sim.dt);
        target = getStimulus(targetStim, targetSOpt);

        %%%
        [weights, mse, rnmse, y] = NLT(target.Signal, nwV, nds);
        
        squMSE(j)      = mse;
        squRNMSE(j) = rnmse;
        squResult(j,:)  = y;     
        
    end
    
    NLTres = struct();
    NLTres.Freqs = Freqs;
    NLTres.Amps = Amps;
    NLTres.sinMSE = sinMSE;
    NLTres.sinRNMSE = sinRNMSE;
    NLTres.sinResult = sinResult;
    NLTres.dblMSE = dblMSE;
    NLTres.dblRNMSE = dblRNMSE;
    NLTres.dblResult = dblResult;
    NLTres.phsMSE = phsMSE;
    NLTres.phsRNMSE = phsRNMSE;
    NLTres.phsResult = phsResult;
    NLTres.squMSE = squMSE;
    NLTres.squRNMSE = squRNMSE;
    NLTres.squResult = squResult;


end

