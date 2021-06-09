function sim = importSim(Components, Stimulus, T, num, saveFolder, nameComment, importSwitch, importStateOnly)
%same params struct as used for multiRun
%need to give stimulus and values for all components free paramaters +
%component type
%num is number at end of the filename. -1 input is for no number
%importSwitch = false then only imports network data when saving in hdf file
%format
%
%
%
% Written by Joel Hochstetter
    
    
    if nargin == 6
       importSwitch = true; 
       importStateOnly = false;
    end
    
    if num == -1 
        num = '';
    else
        num = num2str(num);
    end
    
    switch(Components.ComponentType)
        case 'atomicSwitch' 
            swType = 'a'; 
        case 'tunnelSwitchL'
            swType = 'tl';                      
    end
    
    switch Stimulus.BiasType
        case 'DC'
            Stimulus.stimName   = strcat(Stimulus.BiasType,num2str( Stimulus.Amplitude,3),'V');
        case 'AC'
            Stimulus.stimName = strcat(Stimulus.BiasType, num2str( Stimulus.Amplitude,3),'V_f',num2str(Stimulus.Frequency,3),'Hz');
        case 'ACsaw'
            Stimulus.stimName = strcat(Stimulus.BiasType, num2str( Stimulus.Amplitude,3),'V_f',num2str(Stimulus.Frequency,3),'Hz');
        case 'DCandWait'
           Stimulus.stimName   = strcat(Stimulus.BiasType,num2str(Stimulus.AmplitudeOn,3),'V_off',num2str(Stimulus.OffTime,3),'s', ...
               'offV', num2str(Stimulus.AmplitudeOff,3), 'V');
        case 'DCsaw'
            Stimulus.stimName = strcat(Stimulus.BiasType, num2str( Stimulus.Amplitude,3),'V_f',num2str(Stimulus.Frequency,3),'Hz');  
        case 'Ramp'
            Stimulus.stimName = strcat(Stimulus.BiasType, 'max', num2str( Stimulus.AmplitudeMax,3),'V_min',num2str(Stimulus.AmplitudeMin,3),'V');      
        case 'Custom'
            Stimulus.stimName = 'custom';
        case 'Square'
            Stimulus.stimName   = strcat(Stimulus.BiasType,num2str(Stimulus.AmplitudeOn,3),'V_off',num2str(Stimulus.OffTime,3),'s', ...
               'offV', num2str(Stimulus.AmplitudeOff,3), 'V');
    end
    
    filename = strcat(saveFolder, '/',swType,'_T',num2str(T),'_',Stimulus.stimName,'_s', ...
        num2str(Components.setVoltage,3), '_r', num2str(Components.resetVoltage,3),'_c', ...
        num2str(Components.criticalFlux,3), '_m', num2str(Components.maxFlux,3), '_b', ...
        num2str(Components.boost,3),'_p',num2str(Components.penalty,3), nameComment, num, '.mat')

    sim = load(filename);
    sim = sim.sim;
    
    
    if isfield(sim, 'hdfFile') && importSwitch
        sim.swLam =  h5read(sim.hdfFile, '/swLam');
        if ~importStateOnly
            sim.swV   =  h5read(sim.hdfFile, '/swV');    
            sim.swC   =  h5read(sim.hdfFile, '/swC');  
        end
    end
    
    if isfield(sim.Comp, 'onG') && ~isfield(sim.Comp, 'onR')
        sim.Comp.onR = sim.Comp.onG;
    end
 
    if isfield(sim.Comp, 'offG') && ~isfield(sim.Comp, 'offR')
        sim.Comp.offR = sim.Comp.offG;
    end    
    
    sim.filename = filename;
    
end