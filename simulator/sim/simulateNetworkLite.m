function [OutputDynamics, SimulationOptions] = simulateNetworkLite(Connectivity, Components, Signals, SimulationOptions, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sames as simulateNetworkLite but does not save junction data. Optional to
% save filament state and events (as defined by a threshold)
%
% Left the API of snapshots. For later usage of visualize the network.
% ARGUMENTS: 
% Connectivity - The Connectivity information of the network. Most
%                importantly the edge list, which shows the connectivity
%                condition between nanowires, and number of nodes and
%                junctions.
% Components - Structure that contains the component properties. Every 
%              field is a (E+1)x1 vector. The extra component is the
%              tester resistor connected in series to the voltage and to 
%              the network.
% Stimulus - Structure that contains the details of the external stimulus
%            (time axis details, voltage signal).
% SimulationOptions - Structure that contains general simulation details that are indepedent of 
%           the other structures (eg, dt and simulation length);
% varargin - if not empty, contains an array of indidces in which a
%            snapshot of the conductances and voltages in the network is
%            requested. This indices are based on the length of the simulation.
% OUTPUT:
% OutputDynamics -- is a struct with the activity of the network
%                    .networkConductance - the conductance of the network (between the two 
%                     contacts) as a function of time.
%                    .networkCurrent - the overall current from contact (1) to contact (2) as a
%                     function of time.
% Simulationoptions -- same struct as input, with updated field names
% snapshots - a cell array of structs, holding the conductance and voltage 
%             values in the network, at the requested time-stamps.
        
% REQUIRES:
% updateComponentConductance
% updateComponentState
%
% USAGE:
%{
    Connectivity = getConnectivity(Connectivity);
    contact      = [1,2];
    Equations    = getEquations(Connectivity,contact);
    Components   = initializeComponents(Connectivity.NumberOfEdges,Components)
    Stimulus     = getStimulus(Stimulus);
    
    OutputDynamics = runSimulation(Equations, Components, Stimulus);
%}
%
% Authors:
% Ido Marcus
% Paula Sanz-Leon
% Ruomin Zhu
% Joel Hochstetter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %% Initialize:
    compPtr         = ComponentsPtr(Components);        % using this matlab-style pointer to pass the Components structure by reference
    niterations     = SimulationOptions.NumberOfIterations;
    electrodes      = SimulationOptions.electrodes;
    numOfElectrodes = SimulationOptions.numOfElectrodes;
    E               = Connectivity.NumberOfEdges;
    V               = Connectivity.NumberOfNodes;
    edgeList        = Connectivity.EdgeList.';
    RHS             = zeros(V+numOfElectrodes,1); % the first E entries in the RHS vector.
    LHSinit         = zeros(V+numOfElectrodes, V+numOfElectrodes);
    electrodeCurrent   = zeros(niterations, numOfElectrodes);

    if SimulationOptions.saveFilStateOnly == true
        filamentStates = zeros(niterations, E);
    end

    if SimulationOptions.saveEventsOnly == true
        events  = zeros(niterations, 1);
        componentConductance = compPtr.comp.conductance;
    end
    
    %% Solve equation systems for every time step and update:
    for ii = 1 : niterations
        % Show progress:
        %progressBar(ii,niterations);

        onOrOffOld   = compPtr.comp.conductance > 1.1*Components.offConductance(1);      
        
        % Update conductance values:
        updateComponentConductance(compPtr); 
        componentConductance = compPtr.comp.conductance;
        
        onOrOffNew   = componentConductance > 1.1*Components.offConductance(1);     
        if SimulationOptions.saveEventsOnly == true && SimulationOptions.saveFilStateOnly == false
            events(ii) = sum(abs(onOrOffNew - onOrOffOld));
        end

        % Get LHS (matrix) and RHS (vector) of equation:
        Gmat = zeros(V,V);
        for i = 1:E
            Gmat(edgeList(i,1),edgeList(i,2)) = componentConductance(i);
            Gmat(edgeList(i,2),edgeList(i,1)) = componentConductance(i);
        end

        Gmat = diag(sum(Gmat, 1)) - Gmat;

        LHS          = LHSinit;

        LHS(1:V,1:V) = Gmat;

        for i = 1:numOfElectrodes
            this_elec           = electrodes(i);
            LHS(V+i,this_elec)  = 1;
            LHS(this_elec,V+i)  = 1;
            RHS(V+i)            = Signals{i,1}(ii);
        end

        % Solve equation:
        sol = LHS\RHS;

        tempWireV = sol(1:V);
        compPtr.comp.voltage = tempWireV(edgeList(:,1)) - tempWireV(edgeList(:,2));        
        
        % Update element fields:
        [lambda, ~] = updateComponentState(compPtr, SimulationOptions.dt);    % ZK: changed to allow retrieval of local values
        
        if SimulationOptions.saveFilStateOnly == true
            filamentStates(ii,:) = lambda;
        end
        
        electrodeCurrent(ii,:)   = sol(V+1:end); 

    end
    
    % Calculate network conductance and save:
    OutputDynamics.electrodeCurrent   = electrodeCurrent;
    OutputDynamics.wireVoltage        = sol(1:V)';
    
    OutputDynamics.storevoltage       = compPtr.comp.voltage';
    OutputDynamics.storeCon           = compPtr.comp.conductance';
    
    if SimulationOptions.saveFilStateOnly == true
        if SimulationOptions.saveEventsOnly == true
            d = (Components.criticalFlux(1) - abs(filamentStates))*5/Components.criticalFlux(1);
            d(d < 0.0) = 0.0;
            switchC = tunnelSwitchL(d, 0.81, 0.17, Components.offConductance(1), Components.onConductance(1));
            dG = diff(switchC); dG  = [dG; zeros(1,size(switchC, 2))];
            dGG = abs(dG./switchC)/SimulationOptions.dt;
            events = sum(thresholdCrossingPeaks(dGG, SimulationOptions.EventThreshold),2);
            OutputDynamics.lambda =  compPtr.comp.filamentState';
        else
           OutputDynamics.lambda = filamentStates;
        end
    else
        OutputDynamics.lambda             =  compPtr.comp.filamentState';
    end    
   
    if SimulationOptions.saveEventsOnly == true
        events(1) = 0;
        OutputDynamics.events  = events;
    end

    % Calculate network conductance and save:
    OutputDynamics.networkCurrent    =  electrodeCurrent(:, 2:end);
    OutputDynamics.networkConductance = abs(OutputDynamics.networkCurrent(:,end) ./ Signals{1});


    
end