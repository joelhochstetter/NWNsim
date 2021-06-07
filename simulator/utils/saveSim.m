function [filename] = saveSim(Stimulus,SimulationOptions,Output,Components, Connectivity, filename, runID)
    %saves dataset to a file as a struct
    %assumes that component parameters are the same at all switches

    % save name is generated in genSaveName.m
    
    sim.Stim = Stimulus;
    %sim.Stim = rmfield(sim.Stim,'TimeAxis'); %does not same timevector as this would be double stored
    sim.netC = Output.networkConductance;
    
    if SimulationOptions.megaLiteSave
        sim.Stim = rmfield(sim.Stim,'TimeAxis'); %does not same timevector as this would be double stored
    else
        sim.netI = Output.networkCurrent;
    end
     
    switch(Components.ComponentType)
        case 'atomicSwitch' 
            swType = 'a'; 
        case 'memristor'
            swType = 'm';
        case 'tunnelSwitch'
            swType = 't';
        case 'quantCSwitch'
            swType = 'q';
        case 'tunnelSwitch2'
            swType = 't2';    
        case 'tunnelSwitchL'
            swType = 'tl';            
        case 'linearSwitch'
            swType = 'l';            
        case 'brownModel'
            swType = 'b';                 
    end
    
    sim.T = SimulationOptions.T;
    sim.dt = SimulationOptions.dt;
    sim.seed = SimulationOptions.seed;
    if strcmp(Connectivity.WhichMatrix, 'nanoWires')
        sim.ConnectFile = Connectivity.filename;
    else
        sim.adjMat = Connectivity.weights;
    end
    
    %Save maximal Lyapunov exponent
    if isfield(Output, 'LyapunovMax')
        sim.LyapunovMax = Output.LyapunovMax;
    end
        
    if numel(Connectivity.NewNodes) > 0
        sim.NewNodes = Connectivity.NewNodes;
        sim.NewEdges = Connectivity.NewEdges;        
    end
    
    sim.ContactNodes = SimulationOptions.ContactNodes;
    sim.Comp.setV = Components.setVoltage(1);
    sim.Comp.resetV = Components.resetVoltage(1);
    sim.Comp.pen = Components.penalty;
    sim.Comp.boost = Components.boost;
    sim.Comp.critFlux = Components.criticalFlux(1);
    sim.Comp.maxFlux = Components.maxFlux(1);
    sim.Comp.onG = Components.onConductance(1);
    sim.Comp.offG = Components.offConductance(1);
    sim.Comp.swType = swType;
    sim.Comp.nonpolar = Components.nonpolar;
    sim.Comp.stateEquation = Components.stateEquation;
    sim.Comp.noiseType = Components.noiseType;
    sim.Comp.noiseBeta = Components.noiseBeta;
    sim.Comp.noiseLevel = Components.noiseLevel;
    
    %Save extra parameters if need be
    if isfield(SimulationOptions, 'misc')
        sim.misc = SimulationOptions.misc;
    end
    
    
    %this check occurs again in case of duplicities from parallel sim
    %check if the filename exists already and updates the name 
    if exist(strcat(filename,'.mat'), 'file') 
        %check that the file 
        currFile = load(strcat(filename,'.mat'), 'runID');
        if isfield(currFile, 'runID') && runID == currFile.runID
            %save without issues
        else %increment number        
            num = 1;
            while exist(strcat(filename, '_#', num2str(num), '.mat'), 'file') > 0
                %filename = strcat(filename, num2str(num));
                num = num + 1;
            end
            filename = strcat(filename, '_#', num2str(num));
        end
    end 
    
    if isfield(Output, 'EndTime')
    	sim.EndTime = Output.EndTime;
    end
    
    
    if isfield(Output, 'MaxG')
        sim.MaxG    = Output.MaxG;
    end
    
    %Save switches unless specified not to
    if ~isfield(SimulationOptions, 'saveSwitches') || (isfield(SimulationOptions, 'saveSwitches') && SimulationOptions.saveSwitches)
        if isfield(SimulationOptions, 'hdfSave') && SimulationOptions.hdfSave
            sim.hdfFile = strcat(filename, '.h5');
            h5path = strcat(SimulationOptions.saveFolder, '/', sim.hdfFile);
            %https://au.mathworks.com/help/matlab/ref/h5write.html
            %https://au.mathworks.com/help/matlab/ref/h5create.html
            
            if ~exist(h5path, 'file')
                h5create(h5path,'/swLam', size(Output.lambda)) 
                h5create(h5path,'/swV', size(Output.storevoltage))
                h5create(h5path,'/swC', size(Output.storeCon))
            end
            

            h5write(h5path, '/swLam', Output.lambda)
            h5write(h5path, '/swV', Output.storevoltage)        
            h5write(h5path, '/swC', Output.storeCon)            
            
        else
            sim.swLam = Output.lambda;
            sim.swV   = Output.storevoltage;
            sim.swC   = Output.storeCon; 
            if isfield(Output, 'wireVoltage')
                sim.nwV = Output.wireVoltage;
            end
        end
    else
        sim.finalStates = Output.lambda(end,:);
        if SimulationOptions.saveFilStateOnly && ~SimulationOptions.saveEventsOnly
            sim.hdfFile = strcat(filename, '.h5');
            h5path = strcat(SimulationOptions.saveFolder, '/', sim.hdfFile);
            
            if ~exist(h5path, 'file')
                h5create(h5path,'/swLam', size(Output.lambda)) 
            end              
            h5write(h5path, '/swLam', Output.lambda) 
            
        end
    end    
    
    if isfield(Output, 'events')
        sim.events = Output.events;
    end

    sim.filename = filename;
    
    save(strcat(SimulationOptions.saveFolder, '/', filename,'.mat'), 'sim');
    
    
    
end

