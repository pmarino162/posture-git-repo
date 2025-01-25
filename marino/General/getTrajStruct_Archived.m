function [trajStruct] = getTrajStruct(Data,condFields,trajFields)

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
         trajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
         trajStruct.(['all',trajFieldName]) = [];
         trajStruct.(['avg',trajFieldName]) = [];
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
        %Create trajData for storing results of intermediate operations
        trajData = struct('lengths',[],'allTraj',[],'avgTraj',[]);
        trajData = repmat(trajData,1,numTrajFields);
        %Get trajectory lengths
        numCondTrials = size(tempData,2);
        for trajField = 1:numTrajFields
            trajData(trajField).lengths = zeros(1,numCondTrials);
        end
        for trial = 1:numCondTrials
            for trajField = 1:numTrajFields
               trajFieldName = trajFields{trajField};
               traj = tempData(trial).Traj.(trajFieldName);
               length = size(traj,1);
               trajData(trajField).lengths(trial) = length;
            end
        end
        %Get max traj length and numDims for each trajField, preallocate allTraj
        for trajField = 1:numTrajFields
            %Max length
            lengths = trajData(trajField).lengths;
            maxLength = max(lengths);
            %Num dims
            firstValidTrajInd = find(lengths>0,1,'first');
            trajFieldName = trajFields{trajField};
            traj = tempData(firstValidTrajInd).Traj.(trajFieldName);
            numDims = size(traj,2);
            trajData(trajField).allTraj = NaN(maxLength,numDims,numCondTrials);
        end
        %Store all trajectories
        for trial = 1:numCondTrials
            for trajField = 1:numTrajFields
               trajFieldName = trajFields{trajField};
               traj = tempData(trial).Traj.(trajFieldName);
               timestamps = tempData(trial).Traj.([trajFieldName,'Timestamps']);
               length = size(traj,1);
               trajData(trajField).allTraj(1:length,:,trial) = traj;
               trajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
               trajStruct(structInd).(['all',trajFieldName])(trial).traj = traj;
               trajStruct(structInd).(['all',trajFieldName])(trial).timestamps = timestamps;
            end
        end
        %Take condition averages, store
        for trajField = 1:numTrajFields
            allTraj = trajData(trajField).allTraj;
            lengths = trajData(trajField).lengths;
            meanLength = round(mean(lengths(lengths~=0)));
            avgTraj = nanmean(allTraj,3);
            avgTraj = avgTraj(1:meanLength,:);
            trajFieldName = trajFields{trajField};
            trajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
            trajStruct(structInd).(['avg',trajFieldName]) = avgTraj;
        end
        %Update structInd
        structInd = structInd + 1;
    end

end