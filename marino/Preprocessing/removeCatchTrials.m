function [procData] = removeCatchTrials(procData)

%% Identify catch trials
catchTrials = [];
numTrials = size(procData,2);
for trial = 1:numTrials
    trialName = procData(trial).Overview.trialName;
    switch trialName
        %Earl Multi-Posture Grid Reaching 
        case 'GridReachingCatch' 
            catchTrials = [catchTrials,trial];
        case 'Delayed Center Out Catch 20210621'
            catchTrials = [catchTrials,trial];
        %Earl DCO
        case 'HC_CenterOut'
            stateStr = num2str(procData(trial).TrialData.stateTransitions(1,:));
            if contains(stateStr,'3  6')
                catchTrials = [catchTrials,trial];
            end
    end
end

%% Remove catch trials
procData(catchTrials) = [];

end