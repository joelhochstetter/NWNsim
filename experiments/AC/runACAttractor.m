function runACAttractor(Amps, Freqs, dt, T, attractorFolder, connFile, nameComment)
%   Inputs:
%       Amps (vector of doubles): Amplitudes (in volts) on AC triangular signals
%       Freqs (vector of doubles): Frequency (in Hz) on AC triangular signals
%       dt (double): time-step
%       T (double):  Length of simulation.
%               Need to ensure T*Freqs is an integer (integer number of periods
%               Need to manually check convergence of each attractor by
%               plotting G-V (or I-V curves). Otherwise increase T.
%       attractorFolder (string): name of folder to save attractors to
%       connFile (string): name of file containing the network. Defaults to
%              100 nws, 261 jns
%  nameComment (string): name comment. Defaults to no comment.
%
%  Outputs: 
%       Saves simulation files for each attractor
%
% Written by: Joel Hochstetter


    if nargin < 6
        connFile = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
    end
    
    if nargin < 7
        nameComment = '';
    end


    mkdir(attractorFolder)

    %% Initialise simulation paramaters
    params = struct();

    % Set Simulation Options
    params.SimOpt.saveSim         = true;
    params.SimOpt.useParallel     = true; %can set to true to allow parallel processing
    params.SimOpt.hdfSave         = true;  %saves junction parameters to a 'hdf5' file
    params.SimOpt.saveSwitches         = true; %to save memory set this to false and "hdfSave" to false
    params.SimOpt.stopIfDupName = true; %this parameter only runs simulation if the savename is not used.
    params.SimOpt.T                = T; %length of the simulation in seconds
    params.SimOpt.dt               = dt; %time step
    params.SimOpt.nameComment = nameComment;
    
    
    %Set Stimulus
    params.Stim.BiasType     = 'ACsaw'; % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp' \ 'ACsaw'
    params.Stim.Amplitude    = Amps; 
    params.Stim.Frequency    = Freqs; 

    params.SimOpt.saveFolder = attractorFolder;
    
    %Set Components paramaters
    params.Comp.ComponentType  = 'tunnelSwitchL'; %Set switch model
    params.Comp.onConductance   = 7.77e-5;
    params.Comp.offConductance  = 1e-8;
    params.Comp.setVoltage       = 1e-2;
    params.Comp.resetVoltage   = 1e-2;
    params.Comp.criticalFlux   =  0.01;
    params.Comp.maxFlux        = 0.015;
    params.Comp.penalty        =    1;
    params.Comp.boost          =  10;
    
    %Set connect file
    params.Conn.filename = connFile;    
    

    %% Run simulations
    multiRun(params);
    
end