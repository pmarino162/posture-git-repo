function [avgTraj,timestamps,CI95] = getAvgTraj(trajStruct,trajFieldName,binWidth)

%Gets the average trajectory. If trajField is singleBinFR, just average the
%single observation for each trial together.
           
%% Subselect from trajStruct
upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
trajStruct = trajStruct(1).(['all',upperTrajFieldName]);


if strcmpi(trajFieldName,'singleBinFR')
    allTraj = vertcat(trajStruct.traj);
    avgTraj = mean(allTraj);
    stdTraj = std(allTraj);
    numSamples = size(allTraj,1);    
    timestamps = [];
    CI95 = (1.96.*stdTraj)./sqrt(numSamples);
else
    %% Preallocate avgTraj and timestamps
        numCondTrials = size(trajStruct,2);

        %Get max and min time
        startTimes = nan(1,numCondTrials);
        endTimes = nan(1,numCondTrials);
        for trial = 1:numCondTrials
            trialTimestamps = trajStruct(trial).timestamps;
            if ~isempty(trialTimestamps)
                startTimes(trial) = trialTimestamps(1);
                endTimes(trial) = trialTimestamps(end);
            end
        end
        minTime = min(startTimes,[],'omitnan');
        maxTime = max(endTimes,[],'omitnan');

        %Preallocate
        firstValidTrajInd = 1;
        while isempty(trajStruct(firstValidTrajInd).traj)
            firstValidTrajInd = firstValidTrajInd + 1;
        end      
        numDims = size(trajStruct(firstValidTrajInd).traj,2);
        timestamps = minTime:binWidth:maxTime;
        avgTraj = NaN(length(timestamps),numDims,numCondTrials);

    %% Fill average traj
        for trial = 1:numCondTrials
            trialTimestamps = trajStruct(trial).timestamps;
            trialTraj = trajStruct(trial).traj;
            if ~isempty(trialTimestamps)
                trialStartTime = trialTimestamps(1);
                trialEndTime = trialTimestamps(end);
                trialStartInd = find(timestamps==trialStartTime);
                trialEndInd =  find(timestamps==trialEndTime);
                avgTraj(trialStartInd:trialEndInd,:,trial) = trialTraj;
            end
        end

    %% Take average
        stdTraj = std(avgTraj,0,3,'omitnan');
        numSamples = sum(~isnan(avgTraj),3);
        CI95 = 1.96.*stdTraj./sqrt(numSamples);
        avgTraj = mean(avgTraj,3,'omitnan');
    
    %% Snip to median start and end times
        medStartTime = median(startTimes,'omitnan');
        medEndTime = median(endTimes,'omitnan');
        rmInd = timestamps<medStartTime | timestamps>medEndTime;
        timestamps(rmInd) = [];
        avgTraj(rmInd,:) = [];
        CI95(rmInd,:) = [];
end


    
    
end