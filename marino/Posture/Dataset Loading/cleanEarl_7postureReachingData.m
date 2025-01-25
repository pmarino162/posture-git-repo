function [Data] = cleanEarl_7postureReachingData(Data,dataset)


    switch dataset
        case 'E20210706'
            %Exclude Trials when Earl sat up in chair
            rmTrials = [595:598,806];
            %Exclude start-stop-start reaches
            rmTrials = [rmTrials,256,802];
        case 'E20210707'
            %Exclude Trials when Earl sat up in chair
            rmTrials = [655:658,799];
            %Exclude start-stop-start reaches
            %rmTrials = [rmTrials,];
        case 'E20210708'
            %Exclude Trials when Earl sat up in chair
            rmTrials = [867,937];
            %Exclude start-stop-start reaches
            %rmTrials = [rmTrials,];
        case 'E20210709'
            %Exclude Trials when Earl sat up in chair
            rmTrials = [590:598];
            %Exclude start-stop-start reaches
            %rmTrials = [rmTrials,];
        case 'E20210710'
            %Exclude Trials when Earl sat up in chair
            rmTrials = [];
            %Exclude start-stop-start reaches
            %rmTrials = [rmTrials,];
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