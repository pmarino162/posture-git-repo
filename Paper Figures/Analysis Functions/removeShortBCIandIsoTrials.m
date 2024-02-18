function [Data] = removeShortBCIandIsoTrials(Data,dataset)

%For BCI and isometric tasks, remove any trials for which the entire
%analysis window isn't contained in the center-out movement 
        cutoff = 250;    

        switch dataset
            %For individual tasks, exclude short trials
            case {'E20200316','E20200317','E20200318','E20200319',... %BCI
                    'N20171215','N20180221',...
                    'R20201020','R20201021',...
                    'E20210901','E20190729','E20190830',...
                    'E20200116','E20200117','E20200120',... %Iso
                    }
                kinData = [Data.kinData];
                moveTime = [kinData.moveTime];
                Data = Data(moveTime >= cutoff);   
            %For multiple tasks paradigm, exclude short non-reaching trials
            case {'E20200314'}         
                numTrials = size(Data,2);
                isReachTrial = false(1,numTrials);
                moveTime = zeros(1,numTrials);
                for trial = 1:numTrials
                    trialName = Data(trial).trialName;
                    if strcmpi(trialName,'HC_CenterOut_ForceBar_20200314')
                        isReachTrial(trial) = true;
                    else
                        moveTime(trial) = Data(trial).kinData.moveTime;
                    end
                end
                Data = Data(moveTime >= cutoff | isReachTrial);                
        end 

end