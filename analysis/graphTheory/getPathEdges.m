function PathEdges = getPathEdges(PathNodes, EdgeList)
%Inputs:
%   PathNodes is an array of the indices of the nodes of some path
%   EdgeList is a (2xE) list of edges. It can be found by
%   Connectivity.EdgeList

    PathEdges = zeros(numel(PathNodes) - 1,1);
    
    for i = 1:numel(PathEdges)
    %Find node that is lower numerically
        if PathNodes(i) < PathNodes(i+1)
            l = PathNodes(i);
            u = PathNodes(i+1);
        else
            l = PathNodes(i+1);
            u = PathNodes(i);
        end
        %get indices of edges connected from PathNodes(i)
        f1 = find(EdgeList(1,:) == l);
        f2 = find(EdgeList(2,:) == u);
        PathEdges(i) = intersect(f1, f2);
    end

end