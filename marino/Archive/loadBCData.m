   function [rawData] = loadBCData(animal,date)
    %% Preallocate, get year and month
    rawData = [];
    year = date(1:4);
    month = date(5:6);
    
    %% Get folder names within date directory
    folderStruct = dir(fullfile('S:\Animals\',animal,year,month,date));
    folderNames = {folderStruct.name};
    
    %% Create functions for finding relevant folders
    @(x) contains(x,'gradualTraining')
    @(x) contains(x,'centerOut')
    containsGT = cellfun(@(x) contains(x,'gradualTraining'),folderNames);
    containsCO = cellfun(@(x) contains(x,'centerOut'),folderNames);
    
    %% Find appropriate folders, load translated files from them
    numFolders = size(folderStruct,1);
    for folder = 1:numFolders
        if containsGT(folder) | containsCO(folder)
           fileStruct = dir(fullfile('S:\Animals\',animal,year,month,date,folderNames{folder}));
           fileNames = {fileStruct.name};
           @(x) contains(x,'translated')
           @(x) contains(x,'SI_SW')
           containsTranslated = cellfun(@(x) contains(x,'translated'),fileNames);
           containsSI = cellfun(@(x) contains(x,'SI_SW'),fileNames);
           numFiles = size(fileStruct,1);
           for file = 1:numFiles
               if containsTranslated(file) & containsSI(file)
                   tempData = load(fullfile('S:\Animals\',animal,year,month,date,folderNames{folder},fileNames{file}));
                   rawData = [rawData,tempData.Data];
               end
           end
        end
    end
    clear tempData
    
end