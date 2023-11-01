clear; clc; close all

%% Main Loop
    NDatasetList = {'N20190222','N20190226','N20190227','N20190228','N20190307'};
    RDatasetList = {'R20200221','R20200222'};
    
    acquireCenterTimes = [];
    centerHoldTimes = [];
    targetHoldTimes = [];
    reachTimes = [];
    for datasetList =NDatasetList
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        numTrials = size(Data,2);
        stateNames = Data(1).stateData.stateNames;
        % Get distribution of hold times
        kinData = [Data.kinData];
        
        centerHoldTimes = [centerHoldTimes,[kinData.centerHoldTime]];
        targetHoldTimes = [targetHoldTimes,[kinData.targetHoldTime]];
        reachTimes = [reachTimes,[kinData.reachTime]];
        for trial = 1:numTrials
            trialStartID = min(find([cellfun(@(x) strcmpi(x,'Trial Start'),stateNames)]==1));
            centerStartID = min(find([cellfun(@(x) strcmpi(x,'Center Hold'),stateNames)]==1));
            stateTransitions = Data(trial).stateData.stateTransitions;
            trialStartTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==trialStartID))));
            centerStartTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==centerStartID))));
            acquireCenterTimes = [acquireCenterTimes,centerStartTime-trialStartTime];
            
        end
    end
    
    
            centerHoldTimesList = unique(centerHoldTimes);
            centerHoldTimesRange = [min(centerHoldTimes),max(centerHoldTimes)];
            figure
                histogram(centerHoldTimes);
                ylabel('Count')
                xlabel('Center Hold Time (ms)')
                
            targetHoldTimesList = unique(targetHoldTimes);
            targetHoldTimesRange = [min(targetHoldTimes),max(targetHoldTimes)];
            figure
                histogram(targetHoldTimes);
                ylabel('Count')
                xlabel('Target Hold Time (ms)')
                

            figure
                histogram(reachTimes);
                ylabel('Count')
                xlabel('reach Time (ms)')
                
            figure
                histogram(acquireCenterTimes);
                ylabel('Count')
                xlabel('acquire center times (ms)')
   %For N only
   pctTargetHoldCatch = (sum(targetHoldTimes>425)/length(targetHoldTimes))*100;
   pctCenterHoldCatch = (sum(centerHoldTimes>125)/length(centerHoldTimes))*100;
   
   %Workspace centers
   targetData = [Data.targetData];