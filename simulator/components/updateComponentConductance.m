%function conductance = updateComponentConductance(compPtr)
function [switchChange, conductance] = updateComponentConductance(compPtr)
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
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switchChange = true;
    switch compPtr.comp.type
        case 'memristor'
            % Relevant fields: 
                % charge - time integral over current for every component.
                % memristance properties - details of the (... _/-\_/-\ ...) shape.
                  
            % Even function of charge:
            charge = abs(compPtr.comp.charge);
            
            % Periodic function of charge:
            %charge = bsxfun(@mod, charge, Components.period);

            % Create a _/-\_ shaped memristance (each component with its own
            % specific values):
            conductance = (charge <= compPtr.comp.lowThreshold).*compPtr.comp.offConductance;

            conductance = ...
                conductance + (charge > compPtr.comp.lowThreshold & ...
                              charge <= compPtr.comp.highThreshold).* ...
                (  ...
                   compPtr.comp.offConductance + ...
                   ((compPtr.comp.onConductance-compPtr.comp.offConductance)./(compPtr.comp.highThreshold-compPtr.comp.lowThreshold)) .* ...
                   (charge-compPtr.comp.lowThreshold) ...
                );

            conductance = conductance + (charge >  compPtr.comp.highThreshold).*compPtr.comp.onConductance;

        case 'atomicSwitch'
            oldOnOff = compPtr.comp.OnOrOff;                       
            
            compPtr.comp.OnOrOff = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;
            compPtr.comp.OnOrOff(compPtr.comp.identity == 0) = true; 
            
%             if ~sum(abs(oldOnOff - compPtr.comp.OnOrOff))
%                 switchChange = false;
%                 conductance = compPtr.comp.conductance;
%                 return
%             end
            
            % passive elements (resistors) are always considered as "open" switches
            
            conductance = (~compPtr.comp.OnOrOff) .* compPtr.comp.offConductance;
            conductance = conductance + ...
                         ( compPtr.comp.OnOrOff) .* compPtr.comp.onConductance;
       
                     
        case 'quantCSwitch'
            oldOnOff = compPtr.comp.OnOrOff;                       
            
            compPtr.comp.OnOrOff = floor(abs(compPtr.comp.filamentState) / compPtr.comp.criticalFlux(1));
            compPtr.comp.OnOrOff(compPtr.comp.identity == 0) = true; 
%             
%             if ~sum(abs(oldOnOff - compPtr.comp.OnOrOff))
%                 switchChange = false;
%                 conductance = compPtr.comp.conductance;
%                 return
%             end
            
            % passive elements (resistors) are always considered as "open" switches
            
            conductance = (~compPtr.comp.OnOrOff) .* compPtr.comp.offConductance;
            conductance = conductance + ...
                         ( compPtr.comp.OnOrOff) .* compPtr.comp.onConductance;
                   
                     
            %adding tunnel conductance
        case 'tunnelSwitch'
            V = compPtr.comp.voltage;
            phi = compPtr.comp.barrHeight; %2;          
            d = (compPtr.comp.criticalFlux - abs(compPtr.comp.filamentState))*5/compPtr.comp.criticalFlux(1) + 0.4;
            d(d<0.4)=0.4;
            conductance = tunnelSwitch(V,d,phi,0.4,compPtr.comp.offConductance(1));
            
            if max(abs(V)) > 2.5*phi %Checks conditions for Simmons is valid
                disp('Results may be inaccurate');
            end    
            compPtr.comp.OnOrOff = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;

            
        %adding tunnel conductance
        case 'tunnelSwitch2'
            V = compPtr.comp.voltage;
            phi = compPtr.comp.barrHeight; %2;          
            d =(compPtr.comp.criticalFlux - abs(compPtr.comp.filamentState))*5/compPtr.comp.criticalFlux(1);
            d(d < 0.0) = 0.0;
            A = compPtr.comp.filArea; %0.17
            conductance = tunnelSwitch2(V, d, phi, A, compPtr.comp.offConductance(1), compPtr.comp.onConductance(1));
            compPtr.comp.OnOrOff = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;
            
        case 'tunnelSwitchL'        
            d = (compPtr.comp.criticalFlux - abs(compPtr.comp.filamentState))*5/compPtr.comp.criticalFlux(1);
            d(d < 0.0) = 0.0;
            phi = compPtr.comp.barrHeight; %2;          
            A = compPtr.comp.filArea; %0.17            
            conductance = tunnelSwitchL(d, phi, A, compPtr.comp.offConductance(1), compPtr.comp.onConductance(1));
            compPtr.comp.OnOrOff = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;

        case 'linearSwitch'         
            lambda = abs(compPtr.comp.filamentState);
            lambda(lambda >= compPtr.comp.criticalFlux(1)) = compPtr.comp.criticalFlux(lambda >= compPtr.comp.criticalFlux(1));
            conductance = linearSwitch(lambda, compPtr.comp.criticalFlux(1), compPtr.comp.offConductance(1), compPtr.comp.onConductance(1));
            compPtr.comp.OnOrOff = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;                                    
            
            
        case 'hybridSwitch'  
            V = compPtr.comp.voltage(1:size(compPtr.comp.filamentState));
            phi = 0.8;          
            d = (0.1-abs(compPtr.comp.filamentState))*30+0.4;
            d(d<0.4)=0.4;
            conductance = tunnelSwitch(V,d,phi,0.4,compPtr.comp.offConductance(1));
            
            if max(abs(V)) > 2.5*phi %Checks conditions for Simmons is valid
                'Results may be inaccurate'
            end
            
            compPtr.comp.OnOrOff = floor(abs(compPtr.comp.filamentState) / compPtr.comp.criticalFlux(1));
            compPtr.comp.OnOrOff(compPtr.comp.identity == 0) = true; 
            
            onRes = tunnelSwitch(V,0.4,phi,0.4,compPtr.comp.offConductance(1));
            
            conductance = conductance + ...
                         ( compPtr.comp.OnOrOff) .* onRes;
                     
        case 'resistor'
                conductance = zeros(size(compPtr.comp.identity)); 
                % That's a place-holder. If conductance is not initialized, the 
                % next statement which takes care of passive elements creates a
                % row rather than a column vector.
                
        case 'nonlinearres'
            conductance =  compPtr.comp.voltage.^2+1e-7;
                
        case 'brownModel'
            %updated in updateComponentState.m
            conductance = (compPtr.comp.OnOrOff.*(compPtr.comp.onConductance - compPtr.comp.offConductance)) + compPtr.comp.offConductance;
    end    
    
    % Components that are resistors have conductance 'lowConductance':
    % treated as having negligible conductance
    % regardless of anything else:
    conductance(compPtr.comp.identity == 0) = compPtr.comp.lowConductance;
    
    % Modify the input with the updated conductance values:
    compPtr.comp.conductance = conductance;
end