function [trajStruct] = getTrajStruct(Data,condFields,trajFields,trialInclStates,varargin)

% This function takes in the Data struct, the desired trajectories, and the 
% condition labels, and returns trajStruct, which contains all trajectories 
% for each condition and condition-averaged trajectories 
%% Variable arguments
    matchConditions = false;
    assignopts(who,varargin);

%% Create extTrialInclStates
    extTrialInclStates = trialInclStates;
    numTrialTypes = size(extTrialInclStates,2);
    for trialType = 1:numTrialTypes
        numStates = size(extTrialInclStates(trialType).inclStates,2);
        addTimeToBeginning = zeros(1,numStates);
        addTimeToBeginning(1,1) = -1000000;
        addTimeToEnd = zeros(1,numStates);
        addTimeToEnd(1,numStates) = 1000000;
        extTrialInclStates(trialType).addTimeToBeginning = mat2cell(addTimeToBeginning,[1],[ones(1,numStates)]);
        extTrialInclStates(trialType).addTimeToEnd = mat2cell(addTimeToEnd,[1],[ones(1,numStates)]);
    end
    
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
    extTrajStruct = struct();
    for condField = 1:numCondFields
         condFieldName = condFields{condField}{1};
         trajStruct.(condFieldName) = [];
         extTrajStruct.(condFieldName) = [];
    end
    numTrajFields = size(trajFields,2);
    for trajField = 1:numTrajFields
         trajFieldName = trajFields{trajField};
         upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
         trajStruct.(['all',upperTrajFieldName]) = [];
         trajStruct.(['avg',upperTrajFieldName]) = [];
         extTrajStruct.(['all',upperTrajFieldName]) = [];
    end
    trajStruct = repmat(trajStruct,1,numCond);
    extTrajStruct = repmat(extTrajStruct,1,numCond);
    
%% Fill trajStruct
    %Add all trials
    structInd = 1;
    for condInd = 1:numCond
        %Fill in condition information in trajStruct
        cond = condList(condInd,:);
        for condField = 1:numCondFields
            condFieldName = condFields{condField}{1};
            trajStruct(structInd).(condFieldName) = cond(condField);
            extTrajStruct(structInd).(condFieldName) = cond(condField);
        end
        %Get condition data
        tempData = Data(all(condLabels==cond,2));
        numCondTrials = size(tempData,2);        
        %Store all trajectories and timestamps
        for trial = 1:numCondTrials
            for trajField = 1:numTrajFields
               trajFieldName = trajFields{trajField};
               [traj,timestamps] = getStatesTraj(tempData(trial),trialInclStates,trajFieldName);
               [extTraj,extTimestamps] = getStatesTraj(tempData(trial),extTrialInclStates,trajFieldName);
               upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
               trajStruct(structInd).(['all',upperTrajFieldName])(trial).traj = traj;
               trajStruct(structInd).(['all',upperTrajFieldName])(trial).timestamps = timestamps;
               extTrajStruct(structInd).(['all',upperTrajFieldName])(trial).traj = extTraj;
               extTrajStruct(structInd).(['all',upperTrajFieldName])(trial).timestamps = extTimestamps;
            end
        end  
        %Update structInd
        structInd = structInd + 1;
    end
    
    %Find minNumCondTrials=N, keep shortest N trials for each condition
    if matchConditions
        %Find minNumCondTrials = N, using first trajField
        numCondTrials = nan(1,numCond);
        trajField = 1;
        trajFieldName = trajFields{trajField};
        upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
        for condInd = 1:numCond
            numCondTrials(condInd) = size(trajStruct(condInd).(['all',upperTrajFieldName]),2);
        end
        minNumCondTrials = min(numCondTrials);
        %Keep shortest N trials for each condition, using first trajField
        for condInd = 1:numCond
            if numCondTrials(condInd)>minNumCondTrials
               for trajField = 1:numTrajFields
                   trajLengths = nan(1,numCondTrials(condInd));
                   trajFieldName = trajFields{trajField};
                   upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
                   for trial = 1:numCondTrials(condInd)
                      trajLengths(trial) = size(trajStruct(condInd).(['all',upperTrajFieldName])(trial).timestamps,2);
                   end
                   [~,sortInd] = sort(trajLengths);
                   trajStruct(condInd).(['all',upperTrajFieldName]) = trajStruct(condInd).(['all',upperTrajFieldName])(sortInd(1:minNumCondTrials));
               end
            end
        end
    end
        
    %Take condition averages, store
    structInd = 1;
    for condInd = 1:numCond
        for trajField = 1:numTrajFields
            trajFieldName = trajFields{trajField};
            upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
            %Compare output to getAvgTraj
            method = 'binning';
            [avgTraj,timestamps] = getAvgTraj(trajStruct(structInd).(['all',upperTrajFieldName]),extTrajStruct(structInd).(['all',upperTrajFieldName]),method);
            %Save timestamps and average trajectory
            trajStruct(structInd).(['avg',upperTrajFieldName]).traj = avgTraj;
            trajStruct(structInd).(['avg',upperTrajFieldName]).timestamps = timestamps;
        end
        %Update structInd
        structInd = structInd + 1;
    end

end