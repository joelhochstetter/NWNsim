function [filename, alreadyExists] = genSaveName(SimulationOptions, Components, Stimulus)
    alreadyExists  = false;
    
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

    if strcmp(Components.ComponentType,  'brownModel')
        filename =  strcat(swType,'_T',num2str(Stimulus.T),'_',Stimulus.stimName,'_sE', ...
            num2str(Components.setEField(1),3), '_rI', num2str(Components.resetCurrent(1),3), ...
            '_sr', num2str(Components.setRate(1),3), '_rr', num2str(Components.decayRate(1),3), ...
            '_A', num2str(Components.tunRes(1),3), '_bh', num2str(Components.barrHeight(1),3), ...            
            '_dm', num2str(Components.gapDistance(1),3), '_wm', num2str(Components.maxFilWidth(1),3), ...                        
            SimulationOptions.nameComment);
    else
        filename =  strcat(swType,'_T',num2str(Stimulus.T),'_',Stimulus.stimName,'_s', ...
            num2str(Components.setVoltage(1),3), '_r', num2str(Components.resetVoltage(1),3),'_c', ...
            num2str(Components.criticalFlux(1),3), '_m', num2str(Components.maxFlux(1),3), '_b', ...
            num2str(Components.boost,3),'_p',num2str(Components.penalty,3),SimulationOptions.nameComment);
    end
   
    %check if the filename exists already and updates the name 
    if exist(strcat(SimulationOptions.saveFolder, '/', filename,'.mat'), 'file') 
        alreadyExists = true;
        num = 1;
        if SimulationOptions.stopIfDupName == false
            while exist(strcat(SimulationOptions.saveFolder, '/', filename, '_#', num2str(num), '.mat'), 'file') > 0
                %filename = strcat(filename, num2str(num));
                num = num + 1;
            end
            
            filename = strcat(filename, '_#', num2str(num));
        end
    end 
    
end