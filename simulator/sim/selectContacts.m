function SimulationOptions = selectContacts(Connectivity, SimulationOptions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chooses two wires as contacts to which an external voltage bias is
% applied.
%
% ARGUMENTS: 
% Connectivity -  Structure that contains the adjacency matrix and its
%                 properties. The graph described by the adjacency matrix 
%                 is said to have V vertices and E edges.
% SimulationOptions - Structure that contains the relevant options for choosing 
%           contacts. It must contain a field .ContactMode, which should be
%           one of the following:
%           - 'farthest' - choose the two wires whose centers are as far as
%                          possible from one another.
%           - 'specifiedDistance' - choose the two wires for which the
%                                   distance between their centers is as
%                                   close to Options.BiProbeDistance as
%                                   possible.
%           - 'random' - choose two wires in random. If the network type is 
%                       'randAdjMat', that's the only option.
%
% OUTPUT:
% SimulationOptions with additional field ContactNodes - (2 x 1) row vector with the indices of the wires (=nodes) 
%           between which the external voltage is applied.
%
% Contacts are defined as an ordered pair of wire indices, such that the 
% first one is biased by with respect to the second. A later use of the 
% convention that the second one is grounded means that the first one has 
% the voltage of the source while the other is held at voltage zero.
%
% REQUIRES:
% none
%
% Authors:
% Ido Marcus
% Paula Sanz-Leon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Input control:
    if strcmp(Connectivity.WhichMatrix, 'randAdjMat') && ...
            ~(strcmp(SimulationOptions.ContactMode, 'random') || strcmp(SimulationOptions.ContactMode, 'preSet'))
        error('For networks of type "randAdjMat", contactMode "random" or "preSet" must be selected');
    end
    
    % Choose contacts:
    switch SimulationOptions.ContactMode
        case 'farthest'
            % Choosing the two wires that are furthest away from each other:
            [contA,contB] = ind2sub(size(Connectivity.wireDistances),find(Connectivity.wireDistances(:) == max(Connectivity.wireDistances(:)),1));
            SimulationOptions.ContactNodes = [contA,contB];
            SimulationOptions.Source = contA;
            SimulationOptions.Drain = contB;        
            
        case 'specifiedDistance'           
            % Choosing two wires that are about biProbeDistance um apart
            devFromProbe = abs(Connectivity.wireDistances(:) - SimulationOptions.BiProbeDistance);
            [contA,contB] = ind2sub(size(Connectivity.wireDistances), find(devFromProbe == min(devFromProbe),1));
            SimulationOptions.ContactNodes = [contA,contB];
            SimulationOptions.Source = contA;
            SimulationOptions.Drain = contB;     
            
        case 'random'
            % Choosing two contacts in random
            SimulationOptions.ContactNodes = randi(Connectivity.NumberOfNodes,2,1);
            SimulationOptions.Source = SimulationOptions.ContactNodes(1);
            SimulationOptions.Drain = SimulationOptions.ContactNodes(2);   
            
        case 'preSet'
            % Using values that have already been set
            if ~isfield(SimulationOptions, 'ContactNodes')
                error(strcat('AtomicSwitchNetworks:', mfilename,':EmptyField'), ...
              ['Contact nodes have not been defined.']); 
            end

        case 'topoFarthest'
        	pathDist = distances(graph(single(Connectivity.weights)));
            [~, linIdx] = max(pathDist(:));
            %Convert linear index to matrix indices
            [row, col] = ind2sub(size(Connectivity.weights), linIdx);
            SimulationOptions.ContactNodes = [row, col];
            SimulationOptions.Source = row;
            SimulationOptions.Drain = col;   
            
        case 'fixedTopoDistance'
        	pathDist = distances(graph(single(Connectivity.weights)));
            possible = find(pathDist == SimulationOptions.ContactGraphDist);
            if numel(possible) == 0
                error('No contacts of specfied distance'); 
            end
            %then choose the maximum physical distance such that graph
            %distance is picked
            linIdx = find(max(Connectivity.wireDistances(possible)) == Connectivity.wireDistances, 1);
            [row, col] = ind2sub(size(Connectivity.weights), linIdx);
            SimulationOptions.ContactNodes = [row, col];            
            SimulationOptions.Source = row;
            SimulationOptions.Drain = col;               
            
    end
    SimulationOptions.numOfElectrodes = numel(SimulationOptions.ContactNodes);
    
    
end