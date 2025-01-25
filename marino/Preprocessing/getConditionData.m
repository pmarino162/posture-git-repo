function [conditionData] = getConditionData(trialName,trialRawData,trialData)


switch trialName
    %Delayed Center Out w Rewards/Punishment
    case 'CenterOut_20181112' %Prez Choking
        successStateID = min(find([cellfun(@(x) strcmpi(x,'Success'),{trialRawData.Parameters.StateTable.stateName})==1]));
        rewardIntervalName = trialRawData.Parameters.StateTable(successStateID).Interval.name;
        switch rewardIntervalName{1,1}
            case 'Small Reward'
                rewardID = 1;
                rewardName = 'small';
            case 'Medium Reward'
                rewardID = 2;
                rewardName = 'medium';
            case 'Large Reward'
                rewardID = 3;
                rewardName = 'large';
            case 'Jackpot Reward'
                rewardID = 4;
                rewardName = 'jackpot';
        end
        conditionData = struct('rewardName','','rewardID',[]);
        conditionData.rewardName = rewardName;
        conditionData.rewardID = rewardID;
        
    %Earl Multi-Posture Grid Reaching
    case {'GridReaching','GridReachingCatch','DelayedCenterOut20210828'}
        if isfield(trialRawData,'conditionData')
            conditionData = trialRawData.conditionData;
        else
            conditionData = [];
        end
        
    %Earl Isometric Force
    case {'IsometricForce_1D'}
        if isfield(trialRawData,'conditionData')
            conditionData = trialRawData.conditionData;
        else
            conditionData = [];
        end
    
    
%         workspaceCenter = trialData.targetData.workspaceCenter(1:2);
%         
%         if isequal(workspaceCenter,[35,-390])
%                 posture = 'upRight';
%                 postureID = 1;
%         elseif isequal(workspaceCenter,[-65,-390])
%                 posture = 'upLeft';
%                 postureID = 2;
%         elseif isequal(workspaceCenter,[-65,-470])
%                 posture = 'downLeft';
%                 postureID = 3;
%         elseif isequal(workspaceCenter,[35,-470])
%                 posture = 'downRight';
%                 postureID = 4;
%         else
%                 posture = '';
%                 postureID = [];
%         end
%         conditionData = struct('posture','','postureID',[]);
%         conditionData.posture = '';
%         conditionData.postureID = [];
        
    %Earl Multi-Posture BCI
    case {'BCI Center Out'}
        if isfield(trialRawData,'conditionData')
            conditionData = trialRawData.conditionData;
        else
            conditionData = [];
        end
        
    %Earl Multi-Posture Posture Device Active Movmements 
    case {'Posture Device Active Movements'}
        if isfield(trialRawData,'conditionData')
            conditionData = trialRawData.conditionData;
        end
        
    otherwise
        conditionData = [];   
end

end