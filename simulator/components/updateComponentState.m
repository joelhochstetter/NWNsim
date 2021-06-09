function [local_lambda, local_voltage] = updateComponentState(compPtr, dt)
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

    
    compPtr.comp.filamentState = compPtr.comp.filamentState + ...
                                 (abs(compPtr.comp.voltage) > compPtr.comp.setVoltage) .* ...
                                 (abs(compPtr.comp.voltage) - compPtr.comp.setVoltage) .* ...
                                 sign(compPtr.comp.voltage) ...
                                 * dt;


    reset = (compPtr.comp.resetVoltage > abs(compPtr.comp.voltage)) .* ...
                                (compPtr.comp.resetVoltage - abs(compPtr.comp.voltage)) .* ...
                                sign(compPtr.comp.filamentState) * dt * compPtr.comp.boost;
    %Reset to 0
    resetTo0 = reset >= abs(compPtr.comp.filamentState);
    compPtr.comp.filamentState = compPtr.comp.filamentState - reset;
    compPtr.comp.filamentState(resetTo0) = 0;

    compPtr.comp.filamentState (compPtr.comp.filamentState >  compPtr.comp.maxFlux) =  compPtr.comp.maxFlux(compPtr.comp.filamentState >  compPtr.comp.maxFlux);
    compPtr.comp.filamentState (compPtr.comp.filamentState < -compPtr.comp.maxFlux) = -compPtr.comp.maxFlux(compPtr.comp.filamentState < -compPtr.comp.maxFlux);


               

    % local values:
    local_lambda = compPtr.comp.filamentState;
    local_voltage = compPtr.comp.voltage;
    
end