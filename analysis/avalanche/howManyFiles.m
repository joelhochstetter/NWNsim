function numFiles = howManyFiles(importMode, importFolder)
    switch importMode
        case 0 %simulated data
            files = dir(strcat(importFolder, '/*.mat'));
        case 1 %TDMS file, datatype of Adrian's file from labview
            files = dir(strcat(importFolder, '/*.tdms'));
        case 2 %text file - rintaro data format
            files = dir(strcat(importFolder, '/*.txt'));
        case 3 %text file - rintaro data format
            files = dir(strcat(importFolder, '/*.mat'));            
    end 
    numFiles = numel(files);
            
    
end