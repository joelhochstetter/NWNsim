function [local_lambda, local_voltage] = updateComponentState(compPtr, dt)
%function updateComponentState(compPtr, dt) % ZK: replaced above for
%retreival of local values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function updates the state field of the input struct (which is passed
% by reference). (for example, the 'charge' fields for memristors or the
% 'filamentState' field for atomic switches.
%
% ARGUMENTS: 
% compPtr - a pointer to a struct containing the properties and current 
%           state of the electrical components of the network.
% dt - length of current timestep.
%
% OUTPUT:
% none
%
% REQUIRES:
% none
%
% Authors:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    switch compPtr.comp.stateEquation
        case 'thresholdPolar' 
            wasOpen = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;
            compPtr.comp.filamentState = compPtr.comp.filamentState + ...
                                         (abs(compPtr.comp.voltage) > compPtr.comp.setVoltage) .* ...
                                         (abs(compPtr.comp.voltage) - compPtr.comp.setVoltage) .* ...
                                         sign(compPtr.comp.voltage) ...
                                         * dt;


            if compPtr.comp.noiseLevel > 0.0
                compPtr.comp.filamentState = compPtr.comp.filamentState + ...
                    junctionNoise(compPtr.comp.noiseType, compPtr.comp.noiseBeta, ...
                    compPtr.comp.noiseLevel, numel(compPtr.comp.filamentState));
            end

            reset = (compPtr.comp.resetVoltage > abs(compPtr.comp.voltage)) .* ...
                                        (compPtr.comp.resetVoltage - abs(compPtr.comp.voltage)) .* ...
                                        sign(compPtr.comp.filamentState) * dt * compPtr.comp.boost;
           %Reset to 0
           resetTo0 = reset >= abs(compPtr.comp.filamentState);
           compPtr.comp.filamentState = compPtr.comp.filamentState - reset;
           compPtr.comp.filamentState(resetTo0) = 0;

           % Filaments that have just disconnected suffer a blow:
           justClosed = wasOpen & (abs(compPtr.comp.filamentState) < compPtr.comp.criticalFlux);

           compPtr.comp.filamentState(justClosed) = compPtr.comp.filamentState(justClosed) / compPtr.comp.penalty(1);

           compPtr.comp.filamentState (compPtr.comp.filamentState >  compPtr.comp.maxFlux) =  compPtr.comp.maxFlux(compPtr.comp.filamentState >  compPtr.comp.maxFlux);
           compPtr.comp.filamentState (compPtr.comp.filamentState < -compPtr.comp.maxFlux) = -compPtr.comp.maxFlux(compPtr.comp.filamentState < -compPtr.comp.maxFlux);

        case 'thresholdNonpolar'
            %update according to modified memristor equation
            compPtr.comp.filamentState = compPtr.comp.filamentState + ...
                             (abs(compPtr.comp.voltage) > compPtr.comp.setVoltage) .* ...
                             (abs(compPtr.comp.voltage) - compPtr.comp.setVoltage) * dt ...
                             + (compPtr.comp.resetVoltage > abs(compPtr.comp.voltage)) .* ...
                             (abs(compPtr.comp.voltage) - compPtr.comp.resetVoltage) * dt ...
                             * compPtr.comp.boost;

            %If lambda < 0 set to 0
            compPtr.comp.filamentState(compPtr.comp.filamentState < 0) = 0;       

            %If lambda > lambda max set to lambdamax
            compPtr.comp.filamentState (compPtr.comp.filamentState >  compPtr.comp.maxFlux) =  compPtr.comp.maxFlux(compPtr.comp.filamentState >  compPtr.comp.maxFlux);
            
        case 'thresholdUnipolar'
            compPtr.comp.filamentState = compPtr.comp.filamentState + ...
                     (abs(compPtr.comp.voltage) > compPtr.comp.setVoltage) .* ...
                     (abs(compPtr.comp.voltage) - compPtr.comp.setVoltage) * dt ...
                     + (compPtr.comp.resetVoltage > abs(compPtr.comp.voltage)) .* ...
                     (abs(compPtr.comp.voltage) - compPtr.comp.resetVoltage) * dt ...
                     * compPtr.comp.boost;

            %If lambda < 0 set to 0
            compPtr.comp.filamentState(compPtr.comp.filamentState < 0) = 0;       

            %If lambda > lambda max set to lambdamax
            compPtr.comp.filamentState(compPtr.comp.filamentState >  compPtr.comp.maxFlux) =  compPtr.comp.maxFlux(compPtr.comp.filamentState >  compPtr.comp.maxFlux);

            %current = abs(compPtr.comp.voltage).* compPtr.comp.conductance;
            %compPtr.comp.filamentState(current >= fuseCurrent) = compPtr.comp.filamentState(current >= fuseCurrent)/1.1;
            %compPtr.comp.filamentState(abs(compPtr.comp.voltage) >= fuseVoltage) =  compPtr.comp.filamentState(abs(compPtr.comp.voltage) >= fuseVoltage)/1.1;
            power = abs(compPtr.comp.voltage).^2 .* compPtr.comp.conductance;
            compPtr.comp.filamentState(power >= compPtr.comp.fusePower) = compPtr.comp.filamentState(power >= compPtr.comp.fusePower)/compPtr.comp.fuseFactor;
     
        case 'HPnonpolar'
            compPtr.comp.filamentState = compPtr.comp.filamentState + ...
                     (compPtr.comp.conductance./compPtr.comp.onConductance) .* ...
                     abs(compPtr.comp.voltage) * dt - compPtr.comp.resetVoltage ./ ...
                     compPtr.comp.maxFlux .* compPtr.comp.filamentState * dt;
                 
            compPtr.comp.filamentState(compPtr.comp.filamentState < 0) = 0;  
             compPtr.comp.filamentState(compPtr.comp.filamentState >  compPtr.comp.maxFlux) =  compPtr.comp.maxFlux(compPtr.comp.filamentState >  compPtr.comp.maxFlux);

        case 'thresholdQC'
            compPtr.comp.filamentState = compPtr.comp.filamentState + ...
                             (abs(compPtr.comp.voltage) > compPtr.comp.setVoltage) .* ...
                             (abs(compPtr.comp.voltage) - compPtr.comp.setVoltage) .* ...
                             sign(compPtr.comp.voltage) ...
                             * dt;
            compPtr.comp.filamentState = compPtr.comp.filamentState - ...
                             (compPtr.comp.resetVoltage > abs(compPtr.comp.voltage)) .* ...
                             (compPtr.comp.resetVoltage - abs(compPtr.comp.voltage)) .* ...
                             sign(compPtr.comp.filamentState) ...
                             * dt .* compPtr.comp.boost;
                         
                         
        case 'thresholdHybrid'
            wasOpen = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;            
            compPtr.comp.filamentState = compPtr.comp.filamentState - dt * compPtr.comp.boost* compPtr.comp.filamentState;
            compPtr.comp.filamentState = compPtr.comp.filamentState + ...
                                         (abs(compPtr.comp.voltage) > compPtr.comp.setVoltage) .* ...
                                         (abs(compPtr.comp.voltage) - compPtr.comp.setVoltage) .* ...
                                         sign(compPtr.comp.voltage) ...
                                         * dt;

            if compPtr.comp.noiseLevel > 0.0
                compPtr.comp.filamentState = compPtr.comp.filamentState + ...
                    junctionNoise(compPtr.comp.noiseType, compPtr.comp.noiseBeta, ...
                    compPtr.comp.noiseLevel, numel(compPtr.comp.filamentState));
            end

           % Filaments that have just disconnected suffer a blow:
           justClosed = wasOpen & (abs(compPtr.comp.filamentState) < compPtr.comp.criticalFlux);
           compPtr.comp.filamentState(justClosed) = compPtr.comp.filamentState(justClosed) / compPtr.comp.penalty(1);            
            
           compPtr.comp.filamentState (compPtr.comp.filamentState >  compPtr.comp.maxFlux) =  compPtr.comp.maxFlux(compPtr.comp.filamentState >  compPtr.comp.maxFlux);
           compPtr.comp.filamentState (compPtr.comp.filamentState < -compPtr.comp.maxFlux) = -compPtr.comp.maxFlux(compPtr.comp.filamentState < -compPtr.comp.maxFlux);            
            
           
        case 'brownModel'
            onOrOff = compPtr.comp.OnOrOff;
            electricField = abs(compPtr.comp.voltage./compPtr.comp.filamentState);
            current       = abs(compPtr.comp.voltage.*compPtr.comp.conductance);
            
            %update filament length
            compPtr.comp.filamentState(electricField > compPtr.comp.setEField) = ...
                compPtr.comp.filamentState(electricField > compPtr.comp.setEField) - ...
                (electricField(electricField > compPtr.comp.setEField) - compPtr.comp.setEField) ...
                .* compPtr.comp.setRate * dt;
            
            compPtr.comp.filamentState(compPtr.comp.filamentState < 0.0) = 0.0;
            
            %update filament width
            compPtr.comp.filamentWidth(current > compPtr.comp.resetCurrent) = ...
                compPtr.comp.filamentWidth(current > compPtr.comp.resetCurrent) - ...
                (current(current > compPtr.comp.resetCurrent) - compPtr.comp.resetCurrent) ...
                .* compPtr.comp.decayRate * dt; 
            
            compPtr.comp.filamentWidth(compPtr.comp.filamentWidth < 0.0) = 0.0;
            
            compPtr.comp.OnOrOff = (compPtr.comp.filamentState == 0.0);
           
            %switched on:
            switchOn = compPtr.comp.OnOrOff > onOrOff;
            compPtr.comp.filamentWidth(switchOn) = compPtr.comp.maxFilWidth(switchOn);
%             compPtr.comp.conductance(switchOn) = compPtr.comp.onConductance(switchOn);
            
            %switched off:
            switchOff = compPtr.comp.OnOrOff < onOrOff;         
            compPtr.comp.filamentState(switchOff) = compPtr.comp.gapDistance(switchOff);
%             compPtr.comp.conductance(switchOff) = compPtr.comp.offConductance(switchOff);
            
            
    end          
               

    % local values:
    local_lambda = compPtr.comp.filamentState;
    local_voltage = compPtr.comp.voltage;
    
end