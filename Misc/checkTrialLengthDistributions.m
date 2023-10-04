clear; clc; clf; close all    

%% Load data
    dataset = 'E20200316';
    [Data,zScoreParams] = loadData(dataset);
    %Get trajStruct
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',true);
      
    % Get number of timestamps for each trial
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
       numTrials = size(trajStruct(i).allZSmoothFR,2);
       for j = 1:numTrials
            numTimestamps = [numTimestamps,size(trajStruct(i).allZSmoothFR(j).traj,1)];
       end
    end        
    figure
        histogram(numTimestamps)
        xlabel('Number of timestamps')
        ylabel('Number of trials')