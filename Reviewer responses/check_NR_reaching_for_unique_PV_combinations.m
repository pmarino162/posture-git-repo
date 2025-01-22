clear; clc; clf; close all

%% Datasets to include in analysis 
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};         
    isoDatasetList = {'E20200116','E20200117','E20200120'};    

[Data,zScoreParams] = loadData('N20190222');

outputMatrix = extractUniqueCombinations(Data)

%% Local function def
function outputMatrix = extractUniqueCombinations(Data)

    % Extract postureID and visualID
    conditionData = [Data.conditionData];
    postureID = [conditionData.postureID];
    visualID = [conditionData.visualID];

    % Combine postureID and visualID into a single matrix
    combinedIDs = [postureID; visualID];

    % Identify when the combination changes
    changes = [true, any(diff(combinedIDs, 1, 2) ~= 0, 1)];

    % Extract all combinations at change points
    uniquePostureID = postureID(changes);
    uniqueVisualID = visualID(changes);

    % Create the output matrix
    outputMatrix = [uniquePostureID; uniqueVisualID];
end