function [trajStruct] = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,varargin)

% This function takes in the Data struct, the desired trajectories, and the 
% condition labels, and returns trajStruct, which contains all trajectories 
% for each condition and condition-averaged trajectories 

%% Variable arguments
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
               [traj,timestamps] = getStatesTraj(tempData(trial),trialInclStates,trajFieldName,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'rmSort',rmSort);
               upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
               trajStruct(structInd).(['all',upperTrajFieldName])(trial).trialNum = tempData(trial).trialNum;
               trajStruct(structInd).(['all',upperTrajFieldName])(trial).traj = traj;
               trajStruct(structInd).(['all',upperTrajFieldName])(trial).timestamps = timestamps;
            end
        end  
        %Update structInd
        structInd = structInd + 1;
    end
      
    %Take condition averages, store
    structInd = 1;
    for condInd = 1:numCond
        for trajField = 1:numTrajFields
            trajFieldName = trajFields{trajField};
            upperTrajFieldName = [upper(trajFieldName(1)),trajFieldName(2:end)];
            [avgTraj,timestamps,CI95] = getAvgTraj(trajStruct(structInd).(['all',upperTrajFieldName]),binWidth);
            trajStruct(structInd).(['avg',upperTrajFieldName]).traj = avgTraj;
            trajStruct(structInd).(['avg',upperTrajFieldName]).timestamps = timestamps;
            trajStruct(structInd).(['avg',upperTrajFieldName]).CI95 = CI95;
        end
        %Update structInd
        structInd = structInd + 1;
    end

    %Get individual trial averages (across timesteps), store
    if getTrialAverages
        for i = 1:size(trajStruct,2)
            %SmoothFR
            if isfield(trajStruct,'allSmoothFR')
                for j = 1:size(trajStruct(i).allSmoothFR,2)
                   trajStruct(i).allSmoothFR(j).trialAvg = mean(trajStruct(i).allSmoothFR(j).traj); 
                end
            end         
            %ZSmoothFR
            if isfield(trajStruct,'allZSmoothFR')
                for j = 1:size(trajStruct(i).allZSmoothFR,2)
                   trajStruct(i).allZSmoothFR(j).trialAvg = mean(trajStruct(i).allZSmoothFR(j).traj); 
                end
            end
        end
    end
    
    
    
end