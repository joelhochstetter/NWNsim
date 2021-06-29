function DC_vary_seed_by_ensemble(seedIdx, netFolder, netSize, saveFolder, Vstar, T, EventThreshold) 
%{
Produces a simulation file at specified seed, netFolder, netSize and Vstar
Converts to discete events inside simulation file

Note: to run this code in any reasonable length in 
   time the loop over seed should be parallelised

    
    Inputs:
              seedIdx: seed of network to use. should be in range [1, smax]
           netFolder: folder containing network ensemble at fixed density.
                                where density  (in units of nanowires/(um)^2
               netSize: network size as a side-length (L um). Network is LxL
          saveFolder: base folder to save simulations.
                        Actual folder: saveFolder/Vstar<Vstar>/seed<seed>
                    Vstar: Voltage in volts
                          T: duration of each simulation in s
  EventThreshold: Event threshold (r) on junctions. Event is defined as 
                                    |dG/dt|/G > r before returning below r 
                                    Defaults to 1e-3
        
  
    Written by: Joel Hochstetter
%}

    %% event threshold on dG/G for junctions for avalanche analysis
    if nargin < 6
        T = 30;
    end
    
    if nargin < 7
        EventThreshold = 1e-3;
    end
    
    
    %% Extract actual seeds for simulations 
    % for a given electrode pairing some networks at low sizes will have no source-drain paths 
    % in this case (eg at L = 50, or density = 0.06) run 'GetConnectedNWNs'
    % and this will give a mapping between seed and seed of connected nwn

    if exist(strcat2({netFolder, '/conn_lx_', netSize, '.mat'}), 'file')
        load(strcat2({netFolder, '/conn_lx_', netSize, '.mat'}), 'conSeeds')
        actualSeed = conSeeds(seedIdx);
        nets = dir(strcat(netFolder, '/*_seed_', num2str(actualSeed,'%03.f'), '*lx_', num2str(netSize), '*.mat'))';                
    else
        nets = dir(strcat(netFolder, '/*_seed_', num2str(seedIdx - 1,'%03.f'), '*lx_', num2str(netSize), '*.mat'))';
    end
    
    %% Get connectivity file
    connFile = strcat(nets(1).folder, '/', nets(1).name); %this file must exist. but if everything before is done connectly it will

    
    %% Set-up name comment and save-name
    nameComment = strcat2({'_Lx', netSize, '_seed', seedIdx - 1}, '%03.f');
    disp(strcat2({'Seed: ', seedIdx - 1, ', L = ', netSize, ', Vstar = ', Vstar}));
    
   
    saveF1 = strcat(saveFolder, '/Vstar', num2str(Vstar), '/seed', num2str(seedIdx - 1,'%03.f'), '/'); %new save folder name
    mkdir(fullfile(saveF1))
    
    %% Run simulation
    DC_Vsweep(saveF1, Vstar*0.01, T, connFile, 0 , '', -1, true, true, 0.025, -1, false, true, EventThreshold, nameComment)

end