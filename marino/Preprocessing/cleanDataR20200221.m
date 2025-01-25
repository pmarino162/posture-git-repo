function [Data] = cleanDataR20200221(Data)

    rmTrials = [1166];


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