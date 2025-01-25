function [trajStruct] = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,varargin)

% This function takes in the Data struct, the desired trajectories, and the 
% condition labels, and returns trajStruct, which contains all trajectories 
% for each condition and condition-averaged trajectories 

%% Variable arguments
    matchConditions = false;
    getTrialAverages = false;
    zScoreParams = [];
    rmSort = [];
    assignopts(who,varargin);

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
         if strcmp(trajFieldName,'zSmoothFR')
             upperTrajFieldName = upperTrajFieldName(2:end);
         end
         trajStruct.(['all',upperTrajFieldName]) = [];
         trajStruct.(['avg',upperTrajFieldName]) = [];
    end
    trajStruct = repmat(trajStruct,1,numCond);
    
%% Fill trajStruct
    %Add all trials
    structInd = 1;
    for condInd = 1:numCond
        condInd
        %Fill in condition information in trajStruct
        cond = condList(condInd,:);
        for condField = 1:numCondFields
            condFieldName = condFields{condField}{1};
            trajStruct(structInd).(condFieldName) = cond(condField);
        end
        %Get condition data
        tempData = Data(all(condLabels==cond,2));
        numCondTrials = size(tempData,2);        
        %Store all trajectories and timestamps
        for trial = 1:numCondTrials
            %trial
            for trajField = 1:numTrajFields
               trajFieldName = trajFields{trajField};
               [traj,timestamps] = getStatesTraj20220419(tempData(trial),trialInclStates,trajFieldName,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'rmSort',rmSort);
               upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
               if strcmp(trajFieldName,'zSmoothFR')
                    upperTrajFieldName = upperTrajFieldName(2:end);
               end
               trajStruct(structInd).(['all',upperTrajFieldName])(trial).trialNum = tempData(trial).trialNum;
               trajStruct(structInd).(['all',upperTrajFieldName])(trial).traj = traj;
               trajStruct(structInd).(['all',upperTrajFieldName])(trial).timestamps = timestamps;
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
        if strcmp(trajFieldName,'zSmoothFR')
             upperTrajFieldName = upperTrajFieldName(2:end);
        end
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
                   if strcmp(trajFieldName,'zSmoothFR')
                        upperTrajFieldName = upperTrajFieldName(2:end);
                   end
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
        condInd
        for trajField = 1:numTrajFields
            trajFieldName = trajFields{trajField};
            upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
            if strcmp(trajFieldName,'zSmoothFR')
                upperTrajFieldName = upperTrajFieldName(2:end);
            end
            [avgTraj,timestamps,CI95] = getAvgTraj20211210(trajStruct(structInd).(['all',upperTrajFieldName]),binWidth);
            %Save timestamps and average trajectory
            trajStruct(structInd).(['avg',upperTrajFieldName]).traj = avgTraj;
            trajStruct(structInd).(['avg',upperTrajFieldName]).timestamps = timestamps;
            trajStruct(structInd).(['avg',upperTrajFieldName]).CI95 = CI95;
        end
        %Update structInd
        structInd = structInd + 1;
    end

    %Get trial averages, store
    if getTrialAverages
        for i = 1:size(trajStruct,2)
            %Individual trials
            for j = 1:size(trajStruct(i).allSmoothFR,2)
               trajStruct(i).allSmoothFR(j).trialAvg = mean(trajStruct(i).allSmoothFR(j).traj); 
            end
            %Condition average
            trajStruct(i).avgSmoothFR.condAvg = mean(trajStruct(i).avgSmoothFR.traj);
        end
    end
    
end