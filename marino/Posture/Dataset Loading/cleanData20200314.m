function [Data] = cleanData20200314(Data)

%Exclude Trials when Earl started reaching too soon
    rmTrials = [418];

    rmTrialInd = [];
    numTrials = size(Data.Data,2);
    for trial = 1:numTrials
       trialNum = str2num(Data.Data(trial).Overview.trialNumber(6:end)); 
       if ismember(trialNum,rmTrials)
           rmTrialInd = [rmTrialInd,trial];
       end
    end
    Data.Data(rmTrialInd) = [];
    
end