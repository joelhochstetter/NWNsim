function savemat(saveFolder, critResults, events, netC)
    if isstruct(critResults)
        save(strcat(saveFolder,'/critResults.mat'), 'critResults');
    end
    
    if numel(events) > 0 
        save(strcat(saveFolder,'/events.mat'),       'events');
    end
    
    if numel(netC) > 0
        save(strcat(saveFolder,'/netC.mat'),       'netC');
    end
    
end