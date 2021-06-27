function li = calcLyapunov(attractorFolder, Attractor, lyFolder, eps, dt, T, R)
% Run Lyapunov Simulations  and saves to file Ly.mat, Ly.h5
% Inputs
%   Enter the attractor file and folder from which the simulations starts
%       from. These are strings
%   lyFolder (string): is the Folder to save to
%   eps (double): is the size of perturbation to each junction
%   dt (size of time-step); Is the size of time-step
%   T: length of simulation time
%   R: number of junctions to run Lyapunov simulations for
%
%
% Ouputs:
%   li is the junction lyapunov exponents (1xE array)
%   Ly.mat stores junction Lyapunov exponents and mean exponent
%   Ly.h5 stores exponential diverges at each time step and running
%       lyapunov exponents at the end of each period to check for convergence
% idx < 0 loops through all things
%
%
% Note: to run this code in any reasonable length in time the loops over
% attractor and loops over perturbing junction should be parallelised
%
%
%
% Written by Joel Hochsteter



	%% Params
	params = struct();

	% Set Simulation Options
	params.SimOpt.useWorkspace    = true;
	params.SimOpt.saveSim         = true;
	params.SimOpt.takingSnapshots = false;
	params.SimOpt.onlyGraphics    = true; %does not plot anything
	params.SimOpt.compilingMovie  = false;
	params.SimOpt.useParallel     = false;
	params.SimOpt.hdfSave         = true;
	params.SimOpt.lyapunovSim     = false;
    params.SimOpt.NodalAnal       = true;
    params.SimOpt.reserveFilename = true;
    params.SimOpt.stopIfDupName = true;
    
	params.SimOpt.dt               = dt;
	params.SimOpt.nameComment      = '';
    params.SimOpt.saveFolder       = strcat(lyFolder, '/', Attractor, '/');
    mkdir(params.SimOpt.saveFolder)
    
    
    %Initialise paramaters from attractor file
    sim = multiImport(struct('SimOpt', struct('saveFolder', attractorFolder), 'importByName', Attractor, 'importStateOnly', true));
    params.Stim.BiasType = sim{1}.Stim.BiasType;
    params.Stim.Amplitude = sim{1}.Stim.Amplitude;    
    params.Stim.Frequency = sim{1}.Stim.Frequency;   
    params.Comp.setVoltage = sim{1}.Comp.setV;
    params.Comp.resetVoltage = sim{1}.Comp.resetV;
    params.Comp.criticalFlux = sim{1}.Comp.critFlux;
    params.Comp.maxFlux = sim{1}.Comp.maxFlux;
    params.Comp.boost = sim{1}.Comp.boost;
    params.Comp.penalty = sim{1}.Comp.pen;
    params.Comp.onConductance = sim{1}.Comp.onG;
    params.Comp.offConductance = sim{1}.Comp.offG;    
	
	params.SimOpt.ContactMode      = 'preSet';
	params.SimOpt.ContactNodes  = sim{1}.ContactNodes;


    % Get number of junctions in network
    if isfield(sim{1},  'ConnectFile')
        conFile = load(sim{1}.ConnectFile);
        E = conFile.number_of_junctions;
        params.Conn.filename = sim{1}.ConnectFile;
    elseif isfield(sim{1},  'adjMat')
            E = (sum(sum(sim{1}.adjMat))/2); %261    
            params.Conn.WhichMatrix = 'adjMat';
            params.Conn.weights = sim{1}.adjMat;        
    end
    
   

    switch(sim{1}.Comp.swType)
        case 'a' 
            swType =  'atomicSwitch' ; 
        case 'tl'
            swType = 'tunnelSwitchL';            
        case 't2'
            swType = 'tunnelSwitch2';             
    end

    params.Comp.ComponentType = swType; %'tunnelSwitchL' is used in paper
	params.SimOpt.T                = T; %may need to increase if converges more slowly
    
    if isfield(sim{1}, 'swLam')
        initLamda                  = sim{1}.swLam(end,1:E)';
    elseif isfield(sim{1}, 'finalStates')
        initLamda                  = sim{1}.finalStates';
    else
        disp('FAILED');
        return;
    end
        
    params.Comp.filamentState = initLamda;
    clear('swLam', 'sim')

    %% Run unperturbed simulation
    %Check if unperturbed simulation already exists
    files = dir(strcat(params.SimOpt.saveFolder, '/*unperturbed.mat'));
    params.SimOpt.saveFilStateOnly = true;
    if numel(files)
        params.importByName = files(1).name;
         params.importStateOnly = true;
        u = multiImport(params);
        params.SimOpt.unpertFilState = u{1}.swLam(:,1:E);
        params = rmfield(params, 'importByName');
    else
        params.SimOpt.nameComment = strcat('_unperturbed');	
        u = multiRun(params);
        params.SimOpt.unpertFilState = u{1}.Output.lambda(:,1:E); 
    end
    clear('u');        
    params.SimOpt.saveFilStateOnly = false;
    
    
    %% Initialise perturbed simulations
    params.SimOpt.lyapunovSim     = true;
    params.SimOpt.LyEps           = eps;    
    params.misc.LyEps             = eps;
    params.SimOpt.saveSwitches    = false;
    params.SimOpt.hdfSave         = false;
	id                            = eye(E)*eps; %identity matrix
    params.importSwitch           = false; %so we don't worry about single switch details
    params.SimOpt.NodalAnal       = false;
	
    numTSteps = round(params.SimOpt.T/ params.SimOpt.dt);
    
    if nargin < 7 %number of junctions to run simulation for
        R = E;
    end
    gij = zeros(numTSteps, R); %exponential divergence at each time-step for each perturbed junction
    


    %% Run perturbed simulations or import
    
    for i = 1:R
        params.Comp.filamentState = initLamda + id(:,i);
        params.SimOpt.nameComment = strcat('_eps', num2str(eps), '_i', num2str(i,'%03.f'));
        params.misc.perturbID     = i;
        files = dir(strcat(params.SimOpt.saveFolder, '/*', params.SimOpt.nameComment, '.mat'));
        if numel(files) == 0  || isempty(who('-file', strcat(files(1).folder, '/', files(1).name), 'sim')) %checks if contains sim
            if numel(files) > 0 && isempty(who('-file', strcat(files(1).folder, '/', files(1).name), 'sim')) %if savefile doesnt contain sim delete sim
                delete(strcat(files(1).folder, '/', files(1).name))
            end                            
            t =  multiRun(params);
            gij(:,i) = t{1}.Output.LyapunovMax; %get exponential divergence at each time-bin
        else
            sprintf('Import, %d\n', i)                      
            t = multiImport(params);
            gij(:,i) = t{1}.LyapunovMax; %get exponential divergence at each time-bin
        end   
        clear('t')     
    end        

    params.SimOpt = rmfield(params.SimOpt, 'unpertFilState'); %remove unperturbed filament state from params for saving
    
   %% Calculate lyapunov and save
    skipFraction = 1.0;%0.3; %Skip time-steps calculating mean. Speeds up convergence of Lyapunov exponent
    numTStep = round(params.SimOpt.T/ params.SimOpt.dt*(1-skipFraction)); %number of time-steps in simulations
    numT = round(params.SimOpt.T * params.Stim.Frequency*(1-skipFraction)); %number of periods
    lij    = zeros(numT, E); %Lyapunov exponent up to end of each period
    gij1   = gij(round(skipFraction*numT)+1:end,:); %exponential divergence skipping initial time-steps

    %Calculate Lyapunov expoent 
    for j = 1:numT
        lij(j,:) = mean(gij1(1:floor(numTStep/numT)*j,:));
    end
    
    if numel(lij) == 0
        lij = 0;
    end
    
    li = lij(end,:); %Junction Lyapunov exponents
    ml = mean(li); %maximum Lyapunov exponent

    hdfFile = strcat(params.SimOpt.saveFolder, '/LyCalc.h5');

    if ~exist(hdfFile, 'file')
        h5create(hdfFile,'/lij', size(lij)) 
        h5create(hdfFile,'/gij', size(gij))
    end

    h5write(hdfFile, '/lij', lij)
    h5write(hdfFile, '/gij',  gij)               
    
    save(strcat(params.SimOpt.saveFolder, '/LyCalc.mat'), 'hdfFile', 'ml', 'params', 'li');


end