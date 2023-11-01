clear; clc; close all

%% Main Loop
    NDatasetList = {'N20171215','N20180221'};
    RDatasetList = {'R20201020','R20201021'};

    centerHoldTimes = [];
    moveTimes = [];
    targetHoldTimes = [];
    for datasetList = RDatasetList
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        numTrials = size(Data,2);
        tempCenterHoldTimes = nan(1,numTrials);
        tempTargetHoldTimes = nan(1,numTrials);
        tempMoveTimes = nan(1,numTrials);
        stateNames = Data(1).stateData.stateNames;
        switch dataset
            case {'N20171215','N20180221'}  
                cursorFreezeID = min(find([cellfun(@(x) strcmpi(x,'Cursor Freeze'),stateNames)]==1));
                cursorReleaseID = min(find([cellfun(@(x) strcmpi(x,'Cursor Release'),stateNames)]==1));
                targetHoldID = min(find([cellfun(@(x) strcmpi(x,'Target Hold'),stateNames)]==1));
                successID = min(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
            case{'R20201020','R20201021'}
                centerHoldStartID = min(find([cellfun(@(x) strcmpi(x,'HoldA'),stateNames)]==1));
                centerHoldEndID = min(find([cellfun(@(x) strcmpi(x,'React'),stateNames)]==1));
                targetHoldStartID = min(find([cellfun(@(x) strcmpi(x,'Hold'),stateNames)]==1));
                successID = min(find([cellfun(@(x) strcmpi(x,'InterTrial'),stateNames)]==1));  
        end
        for trial = 1:numTrials
            stateTransitions = Data(trial).stateData.stateTransitions;
            centerHoldStartTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==centerHoldStartID))));
            centerHoldEndTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==centerHoldEndID))));
            targetHoldStartTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetHoldStartID))));
            successTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==successID))));
            tempCenterHoldTimes(trial) = centerHoldEndTime-centerHoldStartTime;
            tempMoveTimes(trial) = targetHoldStartTime-centerHoldEndTime;
            tempTargetHoldTimes(trial) = successTime-targetHoldStartTime;
        end
        %Add dataset data to all data
        centerHoldTimes = [centerHoldTimes,tempCenterHoldTimes];
        targetHoldTimes = [targetHoldTimes,tempTargetHoldTimes];
        moveTimes = [moveTimes,tempMoveTimes];

    end
    
   
            figure
                histogram(centerHoldTimes);
                ylabel('Count')
                xlabel('Center Hold Time (ms)')
     
            figure
                histogram(targetHoldTimes);
                ylabel('Count')
                xlabel('Target Hold Time (ms)')     
            
            figure
                histogram(moveTimes);
                ylabel('Count')
                xlabel('Move Time (ms)')
                
                