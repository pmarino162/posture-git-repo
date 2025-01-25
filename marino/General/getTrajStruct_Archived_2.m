function [trajStruct] = getTrajStruct(Data,condFields,trajFields,trialInclStates)

% This function takes in the Data struct, the desired trajectories, and the 
% condition labels, and returns trajStruct, which contains all trajectories 
% for each condition and condition-averaged trajectories 

%% Get list of condition labels
    numTrials = size(Data,2);
    numCondFields = size(condFields,2);
    condLabels = zeros(numTrials,numCondFields);
    for trial = 1:numTrials
        for condField = 1:numCondFields
             condLabels(trial,condField) = getfield(Data(trial),condFields{condField}{2:end});
        end
    end
    condList = unique(condLabels,'rows');
    numCond = size(condList,1);
    
%% Create trajStruct, preallocate
    trajStruct = struct();
    for condField = 1:numCondFields
         condFieldName = condFields{condField}{1};
         trajStruct.(condFieldName) = [];
    end
    numTrajFields = size(trajFields,2);
    for trajField = 1:numTrajFields
         trajFieldName = trajFields{trajField};
         upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
         trajStruct.(['all',upperTrajFieldName]) = [];
         trajStruct.(['avg',upperTrajFieldName]) = [];
    end
    trajStruct = repmat(trajStruct,1,numCond);

%% Fill trajStruct
    structInd = 1;
    for condInd = 1:numCond
        %Fill in condition information in trajStruct
        cond = condList(condInd,:);
        for condField = 1:numCondFields
            condFieldName = condFields{condField}{1};
            trajStruct(structInd).(condFieldName) = cond(condField);
        end
        %Get condition data
        tempData = Data(condLabels(:,1)==cond(1) & condLabels(:,2)==cond(2));
        numCondTrials = size(tempData,2);        
        %Store all trajectories and timestamps
        for trial = 1:numCondTrials
            for trajField = 1:numTrajFields
               trajFieldName = trajFields{trajField};
               traj = tempData(trial).Traj.(trajFieldName);
               timestamps = tempData(trial).Traj.([trajFieldName,'Timestamps']);
               upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
               trajStruct(structInd).(['all',upperTrajFieldName])(trial).traj = traj;
               trajStruct(structInd).(['all',upperTrajFieldName])(trial).timestamps = timestamps;
            end
        end        
        %Take condition averages, store
        for trajField = 1:numTrajFields
            trajFieldName = trajFields{trajField};
            upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
            endTimes = nan(1,numCondTrials);
            allSamplePeriods = [];
            for trial = 1:numCondTrials
                trialTimestamps = trajStruct(structInd).(['all',upperTrajFieldName])(trial).timestamps;
                if ~isempty(trialTimestamps)
                    endTimes(trial) = trialTimestamps(end);
                end
                allSamplePeriods = [allSamplePeriods,diff(trialTimestamps)];
            end
            %Get median end time
            endTime = median(endTimes,'omitnan');
            %Get mode inter-timestamp length 
            samplePeriod = mode(allSamplePeriods);
            %Create timestamp vector for average trajectory 
            numForwardBins = ceil((endTime/samplePeriod)-0.5);
            if addTimeToBeginning ~= 0
                numBackBins = ceil((abs(addTimeToBeginning)/samplePeriod)-1/2);
            else
                numBackBins = 0;
            end
            timestamps = -1*(numBackBins*samplePeriod):samplePeriod:numForwardBins*samplePeriod;
            binEdges = [timestamps-samplePeriod/2,timestamps(end)+samplePeriod/2];
%             binEdges = 0-samplePeriod/2:samplePeriod:endTime+samplePeriod/2;
%             timestamps = 0:samplePeriod:endTime;
            %Preallocate average trajectory
            firstValidTrajInd = find(endTimes>0,1,'first');
            traj = tempData(firstValidTrajInd).Traj.(trajFieldName);
            numDims = size(traj,2);
            avgTraj = zeros(size(timestamps,2),numDims);
            %Step through timestamps and get average trajectory value for
            %each
            for i = 1:size(binEdges,2)-1
                binStart = binEdges(i);
                binEnd = binEdges(i+1);
                numPts = 0;
                for trial = 1:numCondTrials    
                    trialTimestamps = trajStruct(structInd).(['all',upperTrajFieldName])(trial).timestamps;
                    traj = trajStruct(structInd).(['all',upperTrajFieldName])(trial).traj;
                    trajInd = find(trialTimestamps >= binStart & trialTimestamps < binEnd);
                    %If there are multiple samples from a trial within a bin, take
                    %average 
                    if size(trajInd,2) > 1
                        avgTraj(i,:) = mean(traj(trajInd,:)) + avgTraj(i,:);
                        numPts = numPts + 1;
                    elseif size(trajInd,2) == 1
                        avgTraj(i,:) = traj(trajInd,:) + avgTraj(i,:);
                        numPts = numPts + 1;
                    end
                end
                avgTraj(i,:) = avgTraj(i,:)/numPts;
            end
            %Save timestamps and average trajectory
            trajStruct(structInd).(['avg',upperTrajFieldName]).traj = avgTraj;
            trajStruct(structInd).(['avg',upperTrajFieldName]).timestamps = timestamps;
        end
        %Update structInd
        structInd = structInd + 1;
    end

end