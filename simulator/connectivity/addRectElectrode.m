function [Connectivity, ContactNodes, SDpath, src, drn] = addRectElectrode(Connectivity, fraction, xFraction)
%{
    Inputs:
        Connectivity, SimulationOptions
        SDpath: is the number of junctions between source and drain not. If
        treading electrode as resistive switch then add 2 to this number
%}

        if nargin < 3
            xFraction  =1.0;
        end

        src = find((Connectivity.VertexPosition(:,2) >= Connectivity.GridSize(2)*(1-fraction)) & ...
            (Connectivity.VertexPosition(:,1) >= Connectivity.GridSize(1)*(1-xFraction)));
        drn = find((Connectivity.VertexPosition(:,2) <= Connectivity.GridSize(2)*fraction) & ...
            (Connectivity.VertexPosition(:,1) <= Connectivity.GridSize(1)*xFraction));        
        Connectivity.addNodes = {src, drn};
        ContactNodes = double(Connectivity.NumberOfNodes) + [1:2]; %add new nodes
        Connectivity          = getConnectivity(Connectivity);
        sp = graphallshortestpaths(sparse(double(Connectivity.weights)));
        SDpath  = sp(ContactNodes(1), ContactNodes(2)) - 2;
        disp(strcat('Rectangular electrode added with num src = ', num2str(numel(src)), ...
            ', num drn =', num2str(numel(drn)), ', sd dist = ', num2str(SDpath)));
end