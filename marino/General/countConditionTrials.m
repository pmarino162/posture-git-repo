function [trialCounts] = countConditionTrials(Data,condFields,varargin)

% This function takes in the Data struct and the condition labels, 
% and returns trialCounts, which contains the number of trials per
% condition

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
    
%% Create trialCounts, preallocate
    trialCounts = struct();
    for condField = 1:numCondFields
         condFieldName = condFields{condField}{1};
         trialCounts.(condFieldName) = [];
    end

    trialCounts.('numTrials') = [];
    trialCounts.('trialNum') = [];
    trialCounts = repmat(trialCounts,1,numCond);
    
%% Fill trialCounts
    %Add all trials
    structInd = 1;
    for condInd = 1:numCond
        condInd
        %Fill in condition information in trialCounts
        cond = condList(condInd,:);
        for condField = 1:numCondFields
            condFieldName = condFields{condField}{1};
            trialCounts(structInd).(condFieldName) = cond(condField);
        end
        %Get condition data
        tempData = Data(all(condLabels==cond,2));
        numCondTrials = size(tempData,2);    
        trialCounts(structInd).numTrials = numCondTrials;
        trialCounts(structInd).trialNum = [tempData.trialNum];
        %Update structInd
        structInd = structInd + 1;
    end

end