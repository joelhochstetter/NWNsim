function [snapshot] = generateSnapshotFromData(swVolt, swLambda, swCon, critLam, netV, netC, timestamp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates snapshot struct for plotting using
% snapshotToFigure.m
%
% Input the following:
%       swVolt     - voltages of each switch
%       swLambda   - filament states for each list
%       swCon      - conductance for each switch
%       Components - Components pointer
%       netV       - network Voltage
%       netC       - network Conductance
%       timestamp  - timestamp when the snapshot was taken
%
% Output:
%           snapshot struct containing contains voltage, resistance, etc.
%
% Written by Joel Hochstetter
    
    %Creates snapshot struct
    snapshot = struct();
    snapshot.Voltage = swVolt;
    snapshot.OnOrOff = abs(swLambda) > critLam ;
    snapshot.Conductance    = swCon;
    snapshot.filamentState = swLambda;
    snapshot.Timestamp = timestamp;
    snapshot.netC = netC; 
    snapshot.netV = netV;
    snapshot.netI = netV*netC;


end