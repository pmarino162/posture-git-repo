clear; clc; close all

%% Main Loop
    NDatasetList = {'N20190222','N20190226','N20190227','N20190228','N20190307'};
    RDatasetList = {'R20200221','R20200222'};

    centerHoldTimes = [];
    targetHoldTimes = [];
    for datasetList = RDatasetList
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        % Get distribution of hold times
        kinData = [Data.kinData];
        centerHoldTimes = [centerHoldTimes,[kinData.centerHoldTime]];
        targetHoldTimes = [targetHoldTimes,[kinData.targetHoldTime]];
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
                
   %For N only
   pctTargetHoldCatch = (sum(targetHoldTimes>425)/length(targetHoldTimes))*100;
   pctCenterHoldCatch = (sum(centerHoldTimes>125)/length(centerHoldTimes))*100;