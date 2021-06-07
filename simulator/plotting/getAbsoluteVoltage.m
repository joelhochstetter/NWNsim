function absoluteVoltage = getAbsoluteVoltage(swVDrop, connectivity, contact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the absolute voltage on wires by integrating the voltage drops
% along paths that start from the grounded wire (contact(2)).
%
% ARGUMENTS: 
% swVDrop - an array of the voltage drops across a given switch
%
% connectivity - a struct containing the networks graph-representation,
%                as returned by 'getConnectivity'. Must contain
%                connectivity.weights and connectivity.EdgeList
% contact - indices of the nanowires (vertices) to which the 
%           external voltage is connected. The second one is assumed
%           to be grounded.
%
% OUTPUT:
% absoluteVoltage - the voltage on each nanowire, assuming that nanowire 
%                   number contact(2) is grounded.
%
% REQUIRES:
% none
%
% Authors:
% Ido Marcus
% Joel Hochstetter
%
%
% Edit 8/9/19: Replaced snapshot with swVDrop
% Pass in the voltage instead of snapshot
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % allocate memory
    absoluteVoltage = zeros(length(connectivity.weights),1);
    visited = false(length(connectivity.weights),1); 
        % when visited(i)=true, vertex number i is already in the que and 
        % its absolute voltage is already known.
    
    % start a BFS traversal from nanowire number contact(2):
    que = contact(2);
    visited(contact(2)) = true;
    while ~isempty(que)
        % In every iteration the vertex in the head of the que is dealt
        % with:
        currVertex = que(1);        
        
        % Calculate the absolute voltage of all its yet unvisited adjacent
        % vertices:
        adjecentUnvisitedVertices = find(connectivity.weights(:,currVertex) & ~visited);
        for i = 1 : length(adjecentUnvisitedVertices)
            nextVertex = adjecentUnvisitedVertices(i);
            if currVertex < nextVertex
                currEdge =   connectivity.EdgeList(1,:) == currVertex ...
                           & connectivity.EdgeList(2,:) == nextVertex    ;
                absoluteVoltage(nextVertex) = absoluteVoltage(currVertex) - swVDrop(currEdge);
            else
                currEdge =   connectivity.EdgeList(1,:) == nextVertex ...
                           & connectivity.EdgeList(2,:) == currVertex    ;
                absoluteVoltage(nextVertex) = absoluteVoltage(currVertex) + swVDrop(currEdge);
            end
        end
        
        % Append all unvisited adjacent vertices to the list:
        que = [que ; find(connectivity.weights(:,currVertex) & ~visited)]; %#ok<AGROW>
        
        % Update flags:
        visited(connectivity.weights(:,currVertex)==1) = true;
        
        % Remove the head of the que:
        que = que(2:end);
    end
end