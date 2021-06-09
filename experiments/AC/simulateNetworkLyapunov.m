function [OutputDynamics, SimulationOptions] = simulateNetworkLyapunov(Connectivity, Components, Signals, SimulationOptions, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modifies simulate network to perform lyapunov exponent calculation
%
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
% Authors:
% Joel Hochstetter
% Ido Marcus
% Paula Sanz-Leon
% Ruomin Zhu
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

    if SimulationOptions.saveSwitches    
        wireVoltage        = zeros(niterations, V);
        junctionVoltage    = zeros(niterations, E);
        junctionConductance = zeros(niterations, E);
        junctionFilament   = zeros(niterations, E);
    end
  
    %Calculate the unperturbed orbit
    unpertFilState        = SimulationOptions.unpertFilState';
    
    LyapunovMax        = zeros(niterations,1); %exponential divergence at each time-step
    
    %% Solve equation systems for every time step and update:
    for ii = 1 : niterations
        
        % Update conductance values:
        updateComponentConductance(compPtr); 
        componentConductance = compPtr.comp.conductance;
        
        % Get LHS (matrix) and RHS (vector) of equation:
        Gmat = zeros(V);

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
        
        %calculate junction voltages
        tempWireV = sol(1:V);
        compPtr.comp.voltage = tempWireV(edgeList(:,1)) - tempWireV(edgeList(:,2));
        
        
        % Update element fields:
        updateComponentState(compPtr, SimulationOptions.dt);    
        
        %Calculate state difference
        deltaLam     = compPtr.comp.filamentState - unpertFilState(:,ii);
        normDeltaLam = norm(deltaLam);
        if normDeltaLam == 0 %if perturbation disappears then either perturbation is too small or time-step is too big
            LyapunovMax(ii) = -inf;
            break
        end
        
        %Update trajectory
        compPtr.comp.filamentState =  unpertFilState(:,ii) + SimulationOptions.LyEps/normDeltaLam*deltaLam;
        
        %Calculate exponential divergence for time-step
        LyapunovMax(ii) = log(normDeltaLam/SimulationOptions.LyEps);
        
        %Calculate network current
        electrodeCurrent(ii,:)   = sol(V+1:end);

        if SimulationOptions.saveSwitches
            wireVoltage(ii,:)        = sol(1:V);
            junctionVoltage(ii,:)    = compPtr.comp.voltage;
            junctionConductance(ii,:) = compPtr.comp.conductance;
            junctionFilament(ii,:)   = compPtr.comp.filamentState;
        end
    end
    
    if ~SimulationOptions.saveSwitches
        % Calculate network conductance and save:
        OutputDynamics.electrodeCurrent   = electrodeCurrent;
        OutputDynamics.wireVoltage        = sol(1:V)';

        OutputDynamics.storevoltage       = compPtr.comp.voltage';
        OutputDynamics.storeCon           = compPtr.comp.conductance';
        OutputDynamics.lambda             =  compPtr.comp.filamentState';
    else
        % Calculate network conductance and save:
        OutputDynamics.electrodeCurrent   = electrodeCurrent;
        OutputDynamics.wireVoltage        = wireVoltage;

        OutputDynamics.storevoltage       = junctionVoltage;
        OutputDynamics.storeCon           = junctionConductance;
        OutputDynamics.lambda             = junctionFilament;
    end

    % Calculate network conductance and save:
    OutputDynamics.networkCurrent    = electrodeCurrent(:, 2);
    OutputDynamics.networkConductance = abs(OutputDynamics.networkCurrent ./ Signals{1});

    %Save exponential divergences
    OutputDynamics.LyapunovMax       = LyapunovMax/SimulationOptions.dt; %normalise by step size
    
end