function [Data] = cleanRockyBCIData(Data,dataset)

    rmTrials = [];
    
    %Remove trials that immediately failed
    for trial = 1:numel(Data)
        stateTransitions = Data(trial).stateData.stateTransitions;
        if size(stateTransitions,2) < 2
            rmTrials = [rmTrials,trial];
        end
    end    

    rmTrialInd = [];
    numTrials = size(Data,2);
    for trial = 1:numTrials
       trialNum = Data(trial).trialNum; 
       if ismember(trialNum,rmTrials)
           rmTrialInd = [rmTrialInd,trial];
       end
    end
    Data(rmTrialInd) = [];
    
end