function Connectivity = getConnectivity(Connectivity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a structure with an adjacency matrix and its properties.
%
% ARGUMENTS: 
% Connectivity -- A structure containing the options. It must contain a 
%                 field .WhichMatrix, whose value is a string specifying 
%                 the structure to be loaded. Options:
%                 - 'nanoWires' - loads the output of 
%                                 generate_nanowires_network.py.
%                                 Required fields:
%                                   .filename (no default)
%                 - 'randAdjMat' - A random connected graph with a 
%                                  specified number of nodes and average 
%                                  degree.
%                                  Required fields:
%                                    .NumberOfNodes
%                                    .AverageDegree (average number of
%                                                    edges per vertex) 
%                                                   (no default).
%                 - 'NearestNeighbor' - A 1D ring (vertex (i) is
%                                       connected to vertices (i-1) and
%                                       (i+1)).
%                                       Requiread fields:
%                                         .NumberOfNodes
%                 - 'Random' - Same as randAdjMat, but does not guarentee
%                              connectedness, and does support non-binary
%                              and non-symmetric weights. For now, the rest
%                              of the code does not support non-binary and
%                              non-symmetric weights.
%                              Required fields:
%                                .NumberOfNodes
%                                .WeightType
%                                .Symmetric
%                 - Several test cases, see below.
%                 If not specified, the following default values are used:
%                 - .NumberOFNodes = 42
%                 - .Symmetric = true
%                       Possible values = {true, false}; 
%                 - .WeightType = 'binary' 
%                       Possible values = {'binary', 'weighted'}
%
% OUTPUT: 
% Connectivity -- the data is returned in the same structure that is passed 
%                 in. Fields added for all values of .WhichMatrix:
%                 -  .weights - Adjacency matrix.
%                 -  .EdgeList - A 2XE matrix of vertex indices, where 
%                                each column represents an edge. This 
%                                field follows two conventions:
%                                1. edgeList(1,:) < edgeList(2,:)
%                                2. The index of an edge is defined as the 
%                                   index of the column containing the 
%                                   indices of the two vertices it 
%                                   connects.
%                 -  .NumberOfNodes
%                 -  .NumberOfEdges
%                 -  .dx  
%                 -  .NodeStr - A cell array containing strings for 
%                               labelling each region in the matrix.
%                 -  .VertexPosition - Euclidean coordinates for centre of 
%                                      regions, in micrometers.
%                 -  .wireDistances - Matrix of pairwise Euclidean
%                                     distance.
%
%                 For Fields added only when .WhichMatrix='nanoWires',
%                 see below.
%
% REQUIRES:
% spanforest
%
% USAGE:
%{
    % Specify a random matrix with 100 nodes
    Connectivity.WhichMatrix   = 'Random';
    Connectivity.NumberOfNodes = 100;
clc
    % Or specify a pregenerated nano-wires matrix
    Connectivity.WhichMatrix = 'nanoWires';
    Connectivity.filename    = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
    
    %Load it:
    Connectivity = getConnectivity(Connectivity); 
%}

% Authors:
% Ido Marcus
% Paula Sanz-Leon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ~isfield(Connectivity, 'seed')
        Connectivity.seed = 1;
    end
    rng(Connectivity.seed);

    NewNodes = [];
    NewEdges = [];

% Assign default values:
    if ~isfield(Connectivity,'NumberOfNodes')
        Connectivity.NumberOfNodes = 42;
    end
    
    if ~isfield(Connectivity,'Symmetric')
        Connectivity.Symmetric = true;
    end
    
    if ~isfield(Connectivity,'WeightType')
        Connectivity.WeightType = 'binary';
    end

    if ~isfield(Connectivity,'WhichMatrix')
        Connectivity.WhichMatrix = 'nanoWires';
    end
    
    
    %{
        in order to implement different electrode
        must be a cell array of arrays. Each element of cell array is list
        of nodes connected to by this node.
        e.g. {[1,2,6], 8} would add 2 extra nodes. One is connected to
        [1,2,6] and other is connected to 8
    %}
    if ~isfield(Connectivity,'addNodes')
        Connectivity.addNodes = {};
    end    
    
    
% Add fields that are specific to Connectivity types
    switch Connectivity.WhichMatrix    
    %---------------------------------------------------------------------%  

        case 'nanoWires'
            load(Connectivity.filename, 'length_x', 'length_y', 'number_of_wires', 'adj_matrix', 'wire_distances');
            %load(Connectivity.filename, 'xc', 'yc', 'xi', 'yi', 'xa', 'ya', 'xb', 'yb', 'KCLCoeff', 'KVLCoeff');
            load(Connectivity.filename, 'xc', 'yc', 'xi', 'yi', 'xa', 'ya', 'xb', 'yb', 'this_seed');
            
            Connectivity.GridSize       = [length_x,length_y];
            Connectivity.NumberOfNodes  = number_of_wires;
            Connectivity.weights        = adj_matrix;
            Connectivity.wireDistances  = wire_distances;
            Connectivity.VertexPosition = [xc;yc].';
            Connectivity.WireEnds       = [xa;ya;xb;yb].';
            Connectivity.EdgePosition   = [xi;yi].';
            Connectivity.seed = this_seed;

            if exist('KCLCoeff','var')
               Connectivity.KCLCoeff = KCLCoeff;
            end
            
            if exist('KVLCoeff','var')
               Connectivity.KVLCoeff = KVLCoeff; 
            end

    %---------------------------------------------------------------------%  

        case 'adjMat'
            if ~isfield(Connectivity, 'weights')
                load(Connectivity.filename, 'adj_matrix');
                Connectivity.weights        = adj_matrix;
            end
                           
            Connectivity.NumberOfNodes  = size(Connectivity.weights,1);


    end
    
    %add new nodes
    numNew = numel(Connectivity.addNodes);
    if numNew > 0
        oldNumNodes = double(Connectivity.NumberOfNodes);
        Connectivity.NumberOfNodes = oldNumNodes + numNew;
        adjMat = Connectivity.weights;
        Connectivity.weights = (zeros(Connectivity.NumberOfNodes));
        Connectivity.weights(1:oldNumNodes,1:oldNumNodes) = adjMat; %may need to fix this line for non-sparse case
        
        for i = 1:numNew
            Connectivity.weights(oldNumNodes + i, Connectivity.addNodes{i}) = 1;
            Connectivity.weights(Connectivity.addNodes{i}, oldNumNodes + i) = 1;
            if strcmp(Connectivity.WhichMatrix, 'nanoWires')
            %if nanowire rep then places nanowire at centroid and 
                Connectivity.wireDistances(oldNumNodes + i) = 1;
                Connectivity.VertexPosition(oldNumNodes + i,1:2) = mean(Connectivity.VertexPosition(Connectivity.addNodes{i},1:2), 1);
                Connectivity.WireEnds(oldNumNodes + i,1:4) = [Connectivity.VertexPosition(oldNumNodes + i,1:2), Connectivity.VertexPosition(oldNumNodes + i,1:2)] + [-0.5,0,0.5,0];
%                 AddEdges = size(Connectivity.EdgePosition, 1)  + (1:numel(Connectivity.addNodes{i}));
                %Edge position currently doesn't work
                Connectivity.EdgePosition(end + (1:numel(Connectivity.addNodes{i})),1:2) = Connectivity.VertexPosition(Connectivity.addNodes{i},1:2);
            elseif isfield(Connectivity, 'VertexPosition')
                Connectivity.VertexPosition(oldNumNodes + i,1:2) = mean(Connectivity.VertexPosition(Connectivity.addNodes{i},1:2), 1);
            end
%             if isfield(Connectivity,'EdgeList')
%                 Connectivity.EdgeList(1:2, end + (1:numel(Connectivity.addNodes{i}))) = [Connectivity.addNodes{i}; oldNumNodes + i];
%             end
%             NewEdges = [NewEdges, AddEdges];
        end
        NewNodes = [oldNumNodes + (1:numNew)];
        if isfield(Connectivity,'EdgeList')        
            Connectivity = rmfield(Connectivity, 'EdgeList');
        end
    end
    
    
    % Generate fields that are common to all matrices
%     if ~isfield(Connectivity,'EdgeList')
    [ii, jj] = find(tril(Connectivity.weights)); 
    Connectivity.EdgeList = [jj ii]'; 
        % (2XE matrix, must follow the conventions specified above)
%     end
    
    %get new edges
    if numNew > 0
        for i = NewNodes
            [~,NE] = find(Connectivity.EdgeList == i);
            NewEdges = sort([NewEdges;NE]);
        end
    end
    
    Connectivity.NumberOfEdges = size(Connectivity.EdgeList, 2);
 
    
    if ~isfield(Connectivity,'NodeStr')
    % Generate nodes labels
        for ns = 1:Connectivity.NumberOfNodes
            Connectivity.NodeStr{ns} = num2str(ns);
        end
    end
    
    
    if ~isfield(Connectivity,'VertexPosition')
        % Generate position
        Connectivity.VertexPosition = zeros(Connectivity.NumberOfNodes,2);
    end
    
    
    Connectivity.NewNodes = NewNodes;    
    Connectivity.NewEdges = NewEdges;

    
    %% Check connectivity of graph and save 
    Connectivity.SingleComponent = true;
    if numel(unique(conncomp(graph(Connectivity.weights)))) > 1
        Connectivity.SingleComponent = false;
    end
    
    
end
