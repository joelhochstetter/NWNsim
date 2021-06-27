function GetConnectedNWNs(folder, rectFraction)
%{
    for a given electrode pairing some networks at low sizes will have no source-drain paths 
        in this case (eg at L = 50, or density = 0.06) run 'GetConnectedNWNs'
        and this will give a mapping between seed and seed of connected nwn
    
    Works for rectangular electrodes.

    Inputs:
              folder: folder containing networks
    rectFraction: fraction of networks for each electrode


    Outputs:
        for each network size in folder saves a file listing the seeds of
            connected networks called e.g. 'conn_lx_50.mat' for NWN of side
            length 50

    
    Written by Joel Hochstetter

%}


    %% defaults
    if nargin < 2
        rectFraction = 0.025;
    end


    %% Extract seeds and SD path lengths with rectangular electrodes
    files = dir(strcat(folder, '/*asn*.mat'));
    N = numel(files);

    xFraction = 1.0;
    SDdists = zeros(N,1);
    Lx      = zeros(N,1);
    seed = zeros(N,1);

    parfor i = 1:N
        Connectivity = getConnectivity(struct('filename', strcat(folder, '/', files(i).name)));
        [Connectivity, ~, SDpath, ~, ~] = addRectElectrode(Connectivity, rectFraction, xFraction);
        SDdists(i) = SDpath;
        Lx(i) = round(Connectivity.GridSize(1));
        seed(i) = round(Connectivity.seed);
    end
    
    Lvals = sort(unique(Lx));

    %% Save seeds of networks which have source-drain paths
    for L = Lvals'
        conSeeds = sort(seed((SDdists < inf) & (Lx == L)));
        save(strcat2({folder, '/conn_lx_', L, '.mat'}), 'conSeeds')
    end
end