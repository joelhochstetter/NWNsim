function simAvAnalysis(simFolder, saveBaseFolder, Vstar, Lx, density, binSizes, NumSims, Tend, fitML)
%{
    Once DC simulations are run for avalanche analysis this script extracts
        events and conductance time-series before performing avalanche
        analysis using critAnalysis().

    To use must run avalanche files as in runAvalancheSims()

    Inputs:
        simFolder (character vector): is the base folder to find simulation files.
            Must use file structure hierarchy as in runAvalancheSims.
        saveBaseFolder (character vector): the base folder to save avalanches 'critResults.mat'
            files
        Vstar (vector of doubles): the specified voltages to perform
            criticality analysis on
        Lx (vector of doubles): the specified network sizess to perform
            criticality analysis on
        density (vector of doubles): the specified densities to perform
            criticality analysis on
        binSize (vector of doubles): the specified binSizes to use for
            criticality analysis. Negative bin-sizes correspond to
            multiples of <IEI>, positive bin-sizes are in (ms) or whichever
            time-step was used for simulations.
        NumSims (integer): number of simulations for which each Voltage /
            Lx / density. defaults to 1000
        Tend (double): Length of time simulations were run for. Defaults to
            30s.
       fitML (boolean): Use maximum likelihood fitting. Defaults to true


    Output is a file critResults.mat which contains data on IEIs and
        avalanches

    Written by Joel Hochstetter
%}



    %% set defaults
    if nargin < 6
        binSizes = -1;
    end
    
    if nargin < 7
        NumSims = 1000;
    end
    
    if nargin < 8
        Tend = 30;
    end
    
    if nargin < 9
        fitML = true; %uses maximum like-lihood fitting. 
    end
    
    
    %% set-up parameters
    tsteps = round(Tend/1e-3);

    % Note: to speed up run-time set 'fitML = 0'
    %however then avalanche exponents cannot be 
    % trusted only their distributions
    
    %% Loop over Voltage /  Lx / density
    for v = 1:numel(Vstar)
        for l = 1:numel(Lx)
            for d = 1:numel(density)
                V = Vstar(v);
                D = density(d);
                L = Lx(l);
                saveFolder = strcat(saveBaseFolder, '/density', num2str(D, '%.2f'), '/Vstar', num2str(V), '/Lx', num2str(L));    
                mkdir(saveFolder)
                eventsExist = 0;
                eventFolder = strcat(saveFolder, '/events.mat');

                if exist(eventFolder, 'file') > 0
                    eventsExist = 1;
                end
                
                %stores conductance and events for all simulations in long
                % vector. critAnalysis.m extracts only avalanches from
                % inside each simulation so there is no effect of joining
                netC = zeros(NumSims*tsteps, 1);
                events =  zeros(NumSims*tsteps, 1);

                %loop seed
                for s = 0:(NumSims-1)
                    folderName = strcat(simFolder, '/Density', num2str(D, '%.2f'), 'ChangeSize/Vstar', num2str(V), '/seed', num2str(s, '%03.f'));
                    fname = dir(strcat(folderName, '/tl_*Lx', num2str(L,'%03.f'), '_*.mat'));
                    
                    %checks a single file of this specification is found
                    if numel(fname) > 1
                        disp('BROKEN')
                    end

                    if numel(fname) == 0
                        continue
                    end

                    params = struct();
                    params.importByName = fname(1).name;
                    params.SimOpt.saveFolder = folderName;
                    %if events have already been processed then uses events
                    %otherwise processes events directly
                    if eventsExist == 0
                        h5file = dir(strcat(folderName, '/tl_Lx*', num2str(L,'%03.f'), '_*.h5'));
                        if numel(h5file) > 0
                            params.importStateOnly = true;
                            sim = multiImport(params);
                            ddG = findEventsFromFilState(sim{1}, struct('method', 'dGGSpike', 'dGGratio', 1e-3));
                        else 
                            params.importStateOnly = false;
                            sim = multiImport(params);
                           ddG = sim{1}.events;
                        end
                        events((1+s*tsteps):((1+s)*tsteps)) = ddG;
                        netC((1+s*tsteps):((1+s)*tsteps)) = sim{1}.netC(1:tsteps);
                    else
                        params.importSwitch = false;
                        sim = multiImport(params);
                        netC((1+s*tsteps):((1+s)*tsteps)) = sim{1}.netC(1:tsteps);  
                    end
                end

                if eventsExist == 1
                    events = loadmat(eventFolder);
                end
                
                %ensure no events overhanging between simulations
                events(tsteps*[1:NumSims])      = 0;
                events(tsteps*[1:NumSims] - 1) = 0;
                savemat(saveFolder, 0, events, netC); %save events and netC

                %%
                filename = '';                
                timeVec = 1e-3:1e-3:(Tend*NumSims);

                for j = 1:numel(binSizes)
                    binSize = binSizes(j);
                    disp(strcat2({'i: ', v, ',binsize = ', binSize}))
                    saveFolder1 = strcat(saveFolder, '/bs', num2str(binSize), '/');    
                    mkdir(saveFolder1);
                    critResults = critAnalysis(events, 1e-3, netC, timeVec, V, filename, strcat(saveFolder1), fitML, binSize, tsteps);
                    savemat(saveFolder1, critResults, 0, 0);
                end

            end
        end
    end