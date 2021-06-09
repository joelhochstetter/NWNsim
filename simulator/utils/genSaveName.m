function [filename, alreadyExists] = genSaveName(SimulationOptions, Components, Stimulus)
%{
    Generates default save name for simulation

    Inputs:
        SimulationOptions, Components and Stimulus structrs

    Outputs:
                filename: name of fit
           alreadyExists: a boolean indicating whether or not a name in this file
            and folder already exists

    Written by Joel Hochstetter
%}


    alreadyExists  = false;
    
    switch(Components.ComponentType)
        case 'atomicSwitch' 
            swType = 'a'; 
        case 'tunnelSwitchL'
            swType = 'tl';                        
    end

    filename =  strcat(swType,'_T',num2str(Stimulus.T),'_',Stimulus.stimName,'_s', ...
        num2str(Components.setVoltage(1),3), '_r', num2str(Components.resetVoltage(1),3),'_c', ...
        num2str(Components.criticalFlux(1),3), '_m', num2str(Components.maxFlux(1),3), '_b', ...
        num2str(Components.boost,3),'_p',num2str(Components.penalty,3),SimulationOptions.nameComment);
   
    %check if the filename exists already and updates the name 
    if exist(strcat(SimulationOptions.saveFolder, '/', filename,'.mat'), 'file') 
        alreadyExists = true;
        num = 1;
        if SimulationOptions.stopIfDupName == false
            while exist(strcat(SimulationOptions.saveFolder, '/', filename, '_#', num2str(num), '.mat'), 'file') > 0
                num = num + 1;
            end
            
            filename = strcat(filename, '_#', num2str(num));
        end
    end 
    
end