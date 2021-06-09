function conductance = updateComponentConductance(compPtr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function updates the 'conductance' field of the input struct (which is 
% passed by reference).
%
% ARGUMENTS: 
% compPtr - a pointer to a struct containing the properties and current 
%           state of the electrical components of the network.
%
% OUTPUT:
% conductance - conductances of individual switches
% switchChange - true if switches change and false otherwise
%
% REQUIRES:
% none
%
% Authors:
% Ido Marcus, Joel Hochstetter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    switch compPtr.comp.type
        case 'atomicSwitch'
            oldOnOff = compPtr.comp.OnOrOff;                       
            
            compPtr.comp.OnOrOff = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;
            compPtr.comp.OnOrOff(compPtr.comp.identity == 0) = true; 
            
            
            % passive elements (resistors) are always considered as "open" switches
            
            conductance = (~compPtr.comp.OnOrOff) .* compPtr.comp.offConductance;
            conductance = conductance + ...
                         ( compPtr.comp.OnOrOff) .* compPtr.comp.onConductance;
                   
        case 'tunnelSwitchL'        
            d = (compPtr.comp.criticalFlux - abs(compPtr.comp.filamentState))*5/compPtr.comp.criticalFlux(1);
            d(d < 0.0) = 0.0;
            phi = compPtr.comp.barrHeight; %2;          
            A = compPtr.comp.filArea; %0.17            
            conductance = tunnelSwitchL(d, phi, A, compPtr.comp.offConductance(1), compPtr.comp.onConductance(1));
            compPtr.comp.OnOrOff = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;

                     
        case 'resistor'
                conductance = zeros(size(compPtr.comp.identity)); 
                % That's a place-holder. If conductance is not initialized, the 
                % next statement which takes care of passive elements creates a
                % row rather than a column vector.
                
    end    
    
    % Components that are resistors have conductance 'lowConductance':
    % treated as having negligible conductance
    % regardless of anything else:
    conductance(compPtr.comp.identity == 0) = compPtr.comp.lowConductance;
    
    % Modify the input with the updated conductance values:
    compPtr.comp.conductance = conductance;
end