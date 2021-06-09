function [runs] = multify (params)
%Takes in parameters struct classified into following categories:
%SimOpt: SimulationOptions
%Stim:   Stimulus
%Comp:   Components
%Conn:   Connectivity
%Then maps all possible elements to column vectors
%Names of these are stored as a row vector 
%Combinations of all possible parameter sets are then taken
%Simulations are run for each combination
%Generalise so works for generic Struct names
%This works for numeric parameters and for multiple component files
%To enter numeric parameters enter array of possible numbers 
%... e.g. params.Stim.Amplitude = [0.5:0.25:5];
%To enter Components files enter as either as a char array or cell array
%For char array if '.mat' file then loads this file
%Else tries to find a folder of the name entered and runs for all in file
%For cell array all elements must be char arrays ending in '.mat' 
%
%
% Written by Joel Hochstetter


%combs stores the possible combinations of parameters
%the row number corresponds to the variable number
%the field name for a given row number is stored in flds
%The struct to find each field in is given in fldType
%Each column corresponds to a different set of simulation paramaters

    SimOpt  = struct();
    Stim    = struct();
    Comp    = struct();
    Conn    = struct();
    
    flds    = {};
    fldType = [];
    
    combs = [1]; %initial element is a dummy it is deleted at the end
   %isnumeric == 1 if is a number
   %use to classify 
    %creates list of fields for simulation options
    if isfield(params,'SimOpt') > 0 
        soFlds = fieldnames(params.SimOpt);
        for i = 1:numel(soFlds)
            if isnumeric(params.SimOpt.(soFlds{i})) == 0 || strcmp('ContactNodes', soFlds{i})  || strcmp('unpertFilState', soFlds{i})
                SimOpt.(soFlds{i}) = params.SimOpt.(soFlds{i});
            else
                combs = combvec(combs,params.SimOpt.(soFlds{i}));
                flds{end + 1}    = soFlds{i};
                fldType(end + 1) = 1;
            end            
        end
    end
    
    %creates list of fields for stimulus
    if isfield(params,'Stim') > 0     
        stFlds = fieldnames(params.Stim);
        for i = 1:numel(stFlds)
            if isnumeric(params.Stim.(stFlds{i})) == 0 || strcmp('Signal', stFlds{i})
                Stim.(stFlds{i}) = params.Stim.(stFlds{i});
            else
                combs = combvec(combs,params.Stim.(stFlds{i}));
                flds{end + 1}    = stFlds{i};
                fldType(end + 1) = 2;
            end            
        end
    end
    
    
    %creates list of fields for components
    if isfield(params,'Comp') > 0
        cpFlds = fieldnames(params.Comp);
        for i = 1:numel(cpFlds)
            %Initial state vector are specified or we are a lyapunov sim
            if isnumeric(params.Comp.(cpFlds{i})) == 0 || strcmp('filamentState', cpFlds{i})
                Comp.(cpFlds{i}) = params.Comp.(cpFlds{i});
            else
                combs = combvec(combs,params.Comp.(cpFlds{i}));
                flds{end + 1}    = cpFlds{i};
                fldType(end + 1) = 3;
            end            
        end
        
    end
    

    

    %create list of fields for connectivity
    if isfield(params,'Conn') > 0 

 
            
        %Handles multiple file types if enterred as cell array
        if isfield(params.Conn, 'filename')
            %if it is a char array then 
            if iscell(params.Conn.filename)
                 params.Conn.fNamesCells = 1:numel(params.Conn.filename);
            elseif exist(params.Conn.filename, 'dir') 
                %try to find folder and run on all elements in folder                   
                    %connectivity of file
                    params.Conn.filename = dir(strcat(params.Conn.filename, '/*.mat'));
                    if numel(params.Conn.filename) == 0
                        error('Folder specified is empty')
                    end
                    tmp = cell(numel(params.Conn.filename),1);
                    for i = 1:numel(params.Conn.filename)
                        tmp{i} = params.Conn.filename(i).name; 
                    end
                    params.Conn.filename = tmp;
                    params.Conn.fNamesCells = 1:numel(params.Conn.filename);
            end
        end

        cnFlds = fieldnames(params.Conn);
        for i = 1:numel(cnFlds)
            if (isnumeric(params.Conn.(cnFlds{i})) == 0) || strcmp('weights', cnFlds{i})
                Conn.(cnFlds{i}) = params.Conn.(cnFlds{i});
            else
                combs = combvec(combs,params.Conn.(cnFlds{i}));
                flds{end + 1}    = cnFlds{i};
                fldType(end + 1) = 4;
            end            
        end
    end   
        
    combs = combs(2:end,:);
    
    runs  = cell(1,size(combs,2));
    
    for i = 1:size(combs,2)
        runs{i}.SimOpt = SimOpt;
        runs{i}.Stim   = Stim;
        runs{i}.Comp   = Comp;
        runs{i}.Conn   = Conn;
        
        for j = 1:size(combs,1)
            switch fldType(j)
                case 1
                    runs{i}.SimOpt.(flds{j}) = combs(j,i);
                case 2
                    runs{i}.Stim.(flds{j})   = combs(j,i);
                case 3
                    runs{i}.Comp.(flds{j})   = combs(j,i);
                case 4
                    runs{i}.Conn.(flds{j})   = combs(j,i);
                    if strcmp(flds{j}, 'fNamesCells') %sets connectivity file
                        runs{i}.Conn.filename = params.Conn.filename{round(combs(j,i))};
                        runs{i}.Conn = rmfield(runs{i}.Conn,'fNamesCells'); %remove as field                      
                    end
            end 
        end    
    end
    
end