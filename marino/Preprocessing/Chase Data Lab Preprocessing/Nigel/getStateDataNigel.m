function [stateData] = getStateDataNigel(trialName,trialRawData)

%foundCaseFlag == 0 indicates that a case was not found for the trial type.
%foundCaseFlag == 1 indicates that a case was found for the trial type, and 
%   leads to a set of instructions that is common to many trial types. 
%foundCaseFlag == 2 indicates that a case was found for the trial type, but
%   that further instructions will be specific to that trial type.

%% Preallocate stateData
    stateData = struct('stateNames',{},'stateTransitions',[]);
    
%% Match to case
foundCaseFlag = 0;
switch trialName
    %Nigel Posture BC Center Out
    case 'Nigel Posture BC Center Out'
        foundCaseFlag = 1;
        % Assign State Names
        stateNames = {'Start','Cursor Freeze','Cursor Release','Center Exit','Target Hold','Success with Reward','Trial End'};
        stateData(1).stateNames = stateNames;
        % Get State Transitions
        stateTransitions(1,:) = [1:7];
        stateTransitions(2,:) = [0,trialRawData.HoldAStart,trialRawData.HoldAFinish,trialRawData.ReactionFinish,...
            trialRawData.HoldBStart,trialRawData.HoldBFinish,trialRawData.ComputerFinishTime];
        %Round transition times to ms
        stateTransitions(2,:) = round(stateTransitions(2,:));
        %Store
        stateData(1).stateTransitions = stateTransitions;
end


%% Execute Common Instructions; Error message if there isn't a case for the trial type
if foundCaseFlag == 1
    
elseif foundCaseFlag == 0
    fprintf('Didn''t find trial case for states\n')
    stateData = [];
end


end