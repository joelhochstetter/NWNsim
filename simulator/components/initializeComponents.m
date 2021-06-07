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


    % making new and old scripts compatible
    if isfield(Components, 'onResistance')
        disp('Warning: onResistance provided')
        if ~isfield(Components, 'onConductance')
            Components.onConductance = Components.onResistance;
        end
        disp(strcat2({'Using onConductance = ', Components.onConductance}));        
    end

    if isfield(Components, 'offResistance')
        disp('Warning: offResistance provided')
        if ~isfield(Components, 'offConductance')
            Components.offConductance = Components.offResistance;
        end
        disp(strcat2({'Using onConductance = ', Components.offConductance}));        
    end

    if nargin == 2
        NodalAnal = false;
    end

    %initialises to default values if none are given 
    default.onConductance  = 7.77e-5;
    default.offConductance = 1e-8;
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
    %default.OnOrOff       = 0.0;
    default.setVoltage    = 0.3;%1e-2; %0.3
    default.resetVoltage  = 1e-2;%1e-3; %0.01
    default.criticalFlux  = 1e-4;%1e-1; %1e-4
    default.maxFlux       = 0.15;%0.15;
    default.barrHeight    = 0.81; %potential barrier height for tunnelling in V
    default.filArea       = 0.17; %area of filament tip 
    default.penalty       = 1;%10;
    default.boost         = 10;%10;
    default.nonpolar      = false;
    default.unipolar      = false;
    default.hybrid         = false;
    default.fusePower =  2.5e-5; %for uni-polar switch only
    default.fuseFactor = 1.1;
    
    
    
    
    
    if strcmp(Components.ComponentType, 'brownModel')
        %parameters for Brown model
        default.setRate       = 5e-7;
        default.decayRate     = 3e-5;
        default.tunRes        = 1;
        default.gapDistance   = 0.3454;
        default.maxFilWidth   = 1.0;
        default.setEField     = 10;
        default.resetCurrent  = 0.01;
        default.barrHeight    = 100;
        default.onConductance  = 10;     
        default.filamentWidth = default.maxFilWidth;
        default.filamentState  = default.gapDistance;
    end    
    
    
    
    %{    
        Noise can be added to all junctions by the function junctionNoise.m
        Fields are:
            noiseType = 'powerLaw' (1/f^beta noise), 'gaussian'         
            noiseBeta: power law exponent for powe law noise
            noiseLevel: positive number according to size of noise
    %}
    
    
    default.noiseType = 'gaussian';
    default.noiseBeta  = 2;
    default.noiseLevel = 0.0;
    
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
    
    % this describes which equation the code runs by
    %{
        Options are:
        'thresholdNonpolar' : evolution depends only on absolute sign of voltage
        'thresholdUnipolar': nonpolar with off switchings if power exceeds threshold
        'thresholdPolar': original model. Two possible polarities of switch
            are allowed corresponding to filament growth in reverse direction
        'HPnonpolar': an equivalent version to the HP model with rectangular window function
        'thresholdQC': my original bad model for higher conductance quantum
        'brownModel':  Simon Brown's atomic switch model
    
        To implement (as of 26/4):
        'HPbipolar': bipolar HP model. 
        HP model with window functions. e.g. Biolek, Jogeskar, parabolic
    %}
    
    
    if ~isfield(Components, 'stateEquation')
        if Components.nonpolar
            Components.stateEquation = 'thresholdNonpolar';
        elseif Components.unipolar
            Components.stateEquation = 'thresholdUnipolar';
        elseif Components.hybrid
            Components.stateEquation = 'thresholdHybrid';
        elseif strcmp(Components.ComponentType, 'quantCSwitch') || strcmp(Components.ComponentType, 'hybridSwitch')
            Components.stateEquation = 'thresholdQC';
        elseif strcmp(Components.ComponentType, 'brownModel') 
        %Model from Pike et. al 2020 Nano Lettters: doi:10.1021/acs.nanolett.0c01096
        %Described in MyNotes200710
            Components.stateEquation = 'brownModel';
        else
            Components.stateEquation  = 'thresholdPolar';
        end
    end

    
    Components.voltage       = zeros(E,1);             % (Volt)
    Components.conductance    = ones(E,1)*100;             % (Ohm) (memory allocation)
    Components.onConductance  = ones(E,1)*Components.onConductance;   % (Ohm) 1/(12.9 kOhm) = conductance quantum
    Components.offConductance = ones(E,1)*Components.offConductance; %*1e7;   % (Ohm) literature values

    switch Components.ComponentType        
        case 'memristor'
            Components.charge         = zeros(E,1);                                   % (Coulomb)
            % parameters of (... _/-\_/-\ ...) shape:
            Components.lowThreshold   = rand(E,1)*1e-8;                               % (Coulomb) (1V applied across an OFF-state switch will cause it to open up in about 0.1 sec)
            Components.highThreshold  = (1+rand(E,1)*1e1) .*Components.lowThreshold;  % (Coulomb)
            Components.period         = (2+rand(E,1))     .*Components.highThreshold; % (Coulomb) (optional, usually memristance is not a periodic function)
            Components.OnOrOff        = []; % Dummy field only required in atomic switch
        case {'atomicSwitch', 'tunnelSwitch', 'quantCSwitch', 'hybridSwitch', 'tunnelSwitch2', 'tunnelSwitchL', 'linearSwitch'}
            % parameters of filament formation\dissociation:
            Components.setVoltage    = ones(E,1).*Components.setVoltage; %1e-2;    % (Volt) %% sawtooth: 0.3
            Components.resetVoltage  = ones(E,1).*Components.resetVoltage; %1e-3;    % (Volt) %% sawtooth: 0.01
            Components.criticalFlux  = ones(E,1).*Components.criticalFlux; %1e-1;  % (Volt*sec)  %% sawtooth: 1e-4
            Components.maxFlux       = ones(E,1).*Components.maxFlux; %1.5e-1 % (Volt*sec) %% sawtooth: 0.1
            Components.penalty       = Components.penalty; %10
            Components.boost         = Components.boost; %10
            Components.filamentState = ones(E,1) .* Components.filamentState;        % (Volt*sec)
            Components.OnOrOff       = true(E,1); %This gets fixed upon running sim
            
            
        case 'resistor'
            Components.identity      = zeros(E,1);        % 0 for a passive resistor, 1 for an active element
            Components.OnOrOff       = []; % Dummy field only required in atomic swithc
            
        case 'nonlinearres'
           Components.OnOrOff        = [];
           
        case 'brownModel'
            Components.gapDistance   = ones(E,1).*Components.gapDistance;
            Components.filamentState = ones(E,1).*Components.filamentState;
            Components.filamentWidth = ones(E,1).*Components.filamentWidth;
            Components.OnOrOff       = true(E,1); 
            Components.offConductance = Components.tunRes.*exp(Components.gapDistance .* Components.barrHeight);
  
    end
    
    %% add random initial state
    rng(Components.initSeed);
    switch Components.initRandType
        case 'none'
            %do nothing states already initialised 
        case 'binary'
            r = binornd(1, Components.initBinProb, [E,1]);
            Components.filamentState = (1-r)*Components.initRandLower + r*Components.initRandUpper;
        case 'uniform'
            r = rand(E,1);
            Components.filamentState = (1-r)*Components.initRandLower + r*Components.initRandUpper;            
    end
    
    
end