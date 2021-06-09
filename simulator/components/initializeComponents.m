function Components = initializeComponents(E,Components)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializes a structure which holds the characteristics and current state
% of all the electrical elemetns in the network.
%
% ARGUMENTS: 
% E - number of components.
% Components - a structure containing all the options for the components. It
%           must contain a field 'Component Type' which must be one of the
%           following strings:
%           - 'resistor' - passive element
%           - 'memristor' - an element with a charge-dependent conductance
%                           function (memristance) 
%           - 'atomicSwitch' - an element in which switching events are 
%                              driven by voltage.
% 
%
% OUTPUT:
% Components - a struct containing all the properties of the electrical
%              components in the network. {identity, type, voltage, 
%              conductance, onConductance, offConductance} are obligatory 
%              fields, other fields depend on 'componentType'.
%
% REQUIRES:
% none
%
% USAGE:
%{
    Components.ComponentType = 'atomicSwitch'; 
    Components = initializeComponents(Connectivity.NumberOfEdges,Components);
%}
%
% Authors:
% Ido Marcus, Joel Hochstetter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %initialises to default values if none are given 
    default.onConductance  = 7.77e-5;
    default.offConductance = 7.77e-8;
    default.filamentState = 0.0;
    
    %random initial state
    default.initRandType   = 'none'; % {'none', 'binary', 'uniform'}
    default.initRandLower = 0; % lower value for random initial state. For 'binary' and 'uniform'
    default.initRandUpper = 0; % upper value for random initial state. For 'binary' and 'uniform'
    default.initBinProb       = 0.5; %binary probability of being in initRandUpper other than initRandUpper
    default.initSeed           = 0; %random seed for initial seed
    
    default.conductance    = 1e-7; %conductance of passive elements
    default.lowConductance = 0.1; %for passive elements with fixed conductance
    default.passiveRes     = []; %list of passive resistive elements
    default.setVoltage    = 1e-2;
    default.resetVoltage  = 1e-2;
    default.criticalFlux  = 1e-2;
    default.maxFlux       = 1.5e-2;
    default.barrHeight    = 0.81; %potential barrier height for tunnelling in V
    default.filArea       = 0.17; %area of filament tip 
    default.penalty       = 1; %not used in current model
    default.boost         = 10;
    
    fields = fieldnames(default);
    for i = 1:numel(fields)
        if isfield(Components, fields{i}) == 0
            Components.(fields{i}) = default.(fields{i});
        end
    end
      

    Components.identity      = ones(E,1);          % 0 for a passive resistor, 1 for an active element    
    for i = Components.passiveRes
        Components.identity(i) = 0;
    end
    

        % If one wants an element to be active with probability p,
        % Components.identity      = [rand(E,1) <= p ; 0];
        
    Components.type          = Components.ComponentType; % type of active elements ('atomicSwitch' \ 'memristor' \ ...)
    Components.stateEquation  = 'thresholdPolar';
    
    
    Components.voltage       = zeros(E,1);             % (Volt)
    Components.conductance    = ones(E,1)*100;             % (Ohm) (memory allocation)
    Components.onConductance  = ones(E,1)*Components.onConductance;   % (Ohm) 1/(12.9 kOhm) = conductance quantum
    Components.offConductance = ones(E,1)*Components.offConductance; %*1e7;   % (Ohm) literature values

    switch Components.ComponentType        
        case {'atomicSwitch' 'tunnelSwitchL'}
            % parameters of filament formation\dissociation:
            Components.setVoltage    = ones(E,1).*Components.setVoltage; 
            Components.resetVoltage  = ones(E,1).*Components.resetVoltage;
            Components.criticalFlux  = ones(E,1).*Components.criticalFlux;
            Components.maxFlux       = ones(E,1).*Components.maxFlux; 
            Components.penalty       = Components.penalty; 
            Components.boost         = Components.boost; 
            Components.filamentState = ones(E,1) .* Components.filamentState;
            Components.OnOrOff       = true(E,1); %This gets fixed upon running sim  
            
        case 'resistor'
            Components.identity      = zeros(E,1);        % 0 for a passive resistor, 1 for an active element
            Components.OnOrOff       = []; % Dummy field only required in atomic swithc
            
        case 'nonlinearres'
           Components.OnOrOff        = [];
  
    end
    
    
    
end