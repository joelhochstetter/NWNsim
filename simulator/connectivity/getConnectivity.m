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
%                 -  .speed - conduction speed (for time delays in the 
%                             future) in mm/ms or m/s.
%                 -  .dx  
%                 -  .NodeStr - A cell array containing strings for 
%                               labelling each region in the matrix.
%                 -  .VertexPosition - Euclidean coordinates for centre of 
%                                      regions, in micrometers.
%                 -  .wireDistances - Matrix of pairwise Euclidean
%                                     distance.
%                 -  .delay - Matrix of time delays between regions.
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


    %---------------------------------------------------------------------%  
        case 'randAdjMat' 
            temp = rand(Connectivity.NumberOfNodes) < ...
                   Connectivity.AverageDegree/Connectivity.NumberOfNodes;
            temp = tril(temp,-1);
            adjMat = temp + temp';

            % No empty lines:
            for i=find(sum(adjMat)==0)
                next = mod(i,Connectivity.NumberOfNodes)+1;
                adjMat(i,next) = 1;
                adjMat(next,i) = 1;
            end

            % One connected component:
            [F,~] = spanforest(adjMat);
            if length(F) == 1
                disp('Graph happened to be connected.');
            else
                disp('Graph happened NOT to be connected. Connecting...');

                F{end+1} = F {1};
                for i = 1 : length(F)-1
                    prevTree = F{i};
                    currTree = F{i+1};

                    prevVertex = find(sum(prevTree),1);
                    currVertex = find(sum(currTree),1);

                    adjMat(prevVertex,currVertex)=1;
                    adjMat(currVertex,prevVertex)=1;
                end

                [F,~] = spanforest(adjMat);
                if length(F)==1
                    disp('Success!');
                else
                    disp('Failure!');
                end
            end

            Connectivity.weights = adjMat;

    %---------------------------------------------------------------------%  

        case 'NearestNeighbour'
            Connectivity.weights = diag(ones(1,Connectivity.NumberOfNodes-1),1);
            Connectivity.weights = Connectivity.weights + diag(ones(1,Connectivity.NumberOfNodes-1),-1);
            Connectivity.weights(1,end) = 1;
            Connectivity.weights(end,1) = 1;

    %---------------------------------------------------------------------%  

    case 'Lattice'
        if ~isfield(Connectivity, 'sizex')
            Connectivity.sizex = 10;
        end

        if ~isfield(Connectivity, 'sizey')
            Connectivity.sizey = Connectivity.sizex;
        end            

        if ~isfield(Connectivity, 'BondProb') %bond probability
            Connectivity.BondProb = 1;
        end

        if ~isfield(Connectivity, 'RewireProb') %bond probability
            Connectivity.RewireProb = 0;
        end          

        Connectivity.NumberOfNodes = Connectivity.sizex*Connectivity.sizey;
        
        Connectivity.weights = DilutedLattice(Connectivity.sizex, Connectivity.sizey, Connectivity.BondProb, Connectivity.RewireProb, Connectivity.seed);
        
%         nds = 0:(Connectivity.NumberOfNodes - 1);
%         Connectivity.VertexPosition = 1 +  [floor(nds/sizex); mod(nds, sizex)].';              
    %---------------------------------------------------------------------%  

        case 'WattsStrogatz'
            
            if ~isfield(Connectivity, 'beta') 
                Connectivity.beta = 0.0;
            end
            
            if ~isfield(Connectivity, 'EdgesPerNode') 
                Connectivity.EdgesPerNode = 2;
            end
            
              Connectivity.weights = WattsStrogatz(Connectivity.NumberOfNodes, Connectivity.EdgesPerNode, Connectivity.beta, Connectivity.seed);
    
    %---------------------------------------------------------------------%  

        case 'BarabasiAlbert'
            if ~isfield(Connectivity, 'm0') 
                Connectivity.m0 = 2;
            end
            
            if ~isfield(Connectivity, 'm') 
                Connectivity.mm = 2;
            end            
    
            Connectivity.weights = genScaleFree(Connectivity.NumberOfNodes, Connectivity.m0, Connectivity.m, Connectivity.seed);
    %---------------------------------------------------------------------%  
    
        case 'Random'
            if strcmp(Connectivity.WeightType, 'binary')
                Connectivity.weights  = (randi([0, 1], Connectivity.NumberOfNodes, Connectivity.NumberOfNodes)).*(~ eye(Connectivity.NumberOfNodes));
            else
               Connectivity.weights  = (rand(Connectivity.NumberOfNodes)).*(~ eye(Connectivity.NumberOfNodes));
            end

            if Connectivity.Symmetric
                temp = tril(Connectivity.weights,-1);
                Connectivity.weights = temp + temp.';
            end

    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        % 1         3         2         6  
        % o--\/\/\--o--\/\/\--o--\/\/\--o
        % |         4                   |  
        % | _\/\/\_ o _\/\/\-------------
        % |         5                   |
        % | _\/\/\_ o _\/\/\-------------


        case 'TestCase' 
            Connectivity.weights = [0,0,1,1,1,0;
                                    0,0,1,0,0,1;
                                    1,1,0,0,0,0;
                                    1,0,0,0,0,1;
                                    1,0,0,0,0,1;
                                    0,1,0,1,1,0];
            Connectivity.NumberOfNodes = 6;
    %---------------------------------------------------------------------%  

    % This network is intended for debugging purposes
        % o--\/\/\--o--\/\/\--o--\/\/\--o--\/\/\--o--\/\/\--o
        % 1         2         3         4         5         6

        case 'TestLinear'
            Connectivity.NumberOfNodes = 6;
            Connectivity.weights = diag(ones(Connectivity.NumberOfNodes-1, 1), 1) + ...
                                   diag(ones(Connectivity.NumberOfNodes-1, 1), -1); 

    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        % o--\/\/\--o--\/\/\--o
        % 1         2         3 

        case 'TestSeries'
            Connectivity.NumberOfNodes = 3;
            Connectivity.weights = [0 1 0;
                                    1 0 1;
                                    0 1 0];

    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        %               3         4
        %     |--\/\/\--o--\/\/\--o--\/\/\--|
        %  1--o                             o--2
        %     |--\/\/\--o--\/\/\--o--\/\/\--|
        %               5         6

        case 'TestParallel'
            Connectivity.NumberOfNodes = 6;
            Connectivity.weights = [0 0 1 0 1 0;
                                    0 0 0 1 0 1;
                                    1 0 0 1 0 0;
                                    0 1 1 0 0 0;
                                    1 0 0 0 0 1
                                    0 1 0 0 1 0];
                
    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        % o--\/\/\--o--\/\/\--o--\/\/\--o--\/\/\--o--\/\/\--o--\/\/\--o
        % 1         2         3         4         5         6         1

        case 'TestCircular'
            Connectivity.NumberOfNodes = 6;
            Connectivity.weights = diag(ones(Connectivity.NumberOfNodes-1, 1),  1) + ...
                                   diag(ones(Connectivity.NumberOfNodes-1, 1), -1); 
            Connectivity.weights(1, end)   = 1;
            Connectivity.weights(end, 1)   = 1;

    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        % 1         2
        % o--\/\/\--o
        % | \       |
        % \  \      /
        % /   -     \
        % \    -    /
        % /       \ \
        % |        \|
        % o--\/\/\--o
        % 3         4

        case 'TestSquare'
            Connectivity.NumberOfNodes = 4;
            temp = zeros(Connectivity.NumberOfNodes);
            temp(1  ,2:end) = 1;
            temp(2:3,  end) = 1;
            Connectivity.weights = temp + temp';

    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        % 1         2         3
        % o--\/\/\--o--\/\/\--o 
        % |         |         |
        % <         <         <
        % >         >         >
        % <         <         <
        % |4        |5        |6
        % o--\/\/\--o--\/\/\--o 
        % |         |         |
        % <         <         <
        % >         >         >
        % <         <         <
        % |         |         |
        % o--\/\/\--o--\/\/\--o 
        % 7         8         9

        case 'TestRectangular'
            Connectivity.NumberOfNodes = 9;
            temp_vec = [1 1 0 1 1 0 1 1];
            temp = diag(ones(Connectivity.NumberOfNodes-3, 1), 3) + diag(temp_vec, 1);
            Connectivity.weights = temp + temp';
            
        case 'Minimal'
            Connectivity.NumberOfNodes = 2;          
            Connectivity.weights = [0 1; 1 0];

    %---------------------------------------------------------------------%  

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
 
    if ~isfield(Connectivity,'speed')
        Connectivity.speed = 1.0; 
    end
    
    if ~isfield(Connectivity,'dx')
        Connectivity.dx = 1.0; 
    end
    
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
