function [stateData] = getStateDataRocky(trialName,trialTime,trialStateCodes)

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
    %Rocky Posture BC Center Out
    case 'Rocky Posture BC Center Out'
        foundCaseFlag = 1;
        % Assign State Names
        stateNames = {'Center','Hold A','React','Move','Hold','Trial End'};
        stateData(1).stateNames = stateNames;
        % Get State Transitions
        stateTransitions(1,1) = trialStateCodes(1);
        stateTransitions(2,1) = trialTime(1);
        stInd = 2;
        for i = 2:length(trialStateCodes)
            if trialStateCodes(i) ~= trialStateCodes(i-1)
                stateTransitions(1,stInd) = trialStateCodes(i);
                stateTransitions(2,stInd) = trialTime(i);
                stInd = stInd + 1;
            end
        end
        stateData(1).stateTransitions = round(stateTransitions);
end


%% Execute Common Instructions; Error message if there isn't a case for the trial type
if foundCaseFlag == 1
    
elseif foundCaseFlag == 0
    fprintf('Didn''t find trial case for states\n')
    stateData = [];
end


end