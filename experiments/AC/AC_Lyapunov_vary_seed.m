function AC_Lyapunov_vary_seed(seedIdx, netFolder, netSize, saveFolder, Amps, Freqs
%{
Produces a simulation file at specified seed, netFolder, netSize, Amps and
    Freqs. Performs simulations to calculate Lyapunov exponent.
    Simulaneously calculates avalanche statistics.

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
  EventThreshold: Event threshold (r) on junctions. Event is defined as 
                                    |dG/dt|/G > r before returning below r 
                                    Defaults to 1e-3
        
  
    Written by: Joel Hochstetter
%}


    %% Get network names
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

    connFile = nets(1).name; %this file must exist. but if everything before is done connectly it will
    addpath(netFolder);


    %% Run AC attractor
    dt = 5e-4; %set time-step for simulations
    T  = 3000; %Set simulation time: Must run an integer number of periods for all frequencies
    Amps = [5]; %Specify stimulus amplitudes
    Freqs = [0.25, 0.5, 1]; %Specify stimulus frequencies
    attractorFolder = 'S15attractors'; %set save folder for attractor
    for s = 1:500
        runACAttractor(Amps, Freqs, dt, T, attractorFolder, connFile)
    end

    %% Run Lyapunov simulations
    %set parameters for simulations
    eps   = 5e-4; %size of infintesimal perturbation 
    dt    = 5e-4; %size of simulation time-step
    T     = 200;
    R     = 100; %number of junctions to run Lyapunov simulations for


    files = dir(strcat(attractorFolder, '/*.mat')); %get attractor simulation files
    numA = numel(files); %number of attractor files
    lyFolder = 'lyapunov'; %folder to save Lyapunov exponent simulations

    for i = numA
        calcLyapunov(attractorFolder, files(i).name, lyFolder, eps, dt, T, 500)
    end


end