function [Data] = cleanData20210706(Data)

%Exclude Trials when Earl sat up in chair
    rmTrials = [595:598,806];
    
%Exclude start-stop-start reaches
    rmTrials = [rmTrials,256,802];

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