% Excludes trials for which the specificied states last more than 2 
% stdDevs of the mean length

function [Data] = excludeLengths(Data,trialInclStates)
    
    %Create keepStruct
    numTrials = size(Data,2);
    numTrialTypes = size(trialInclStates,2);
    keepStruct = struct('trialType',[],'length',[],'keep',[]);
    keepStruct = repmat(keepStruct,1,numTrials);
    
    %Get type and length of each trial; store in keepStruct
    for trial = 1:numTrials
        %Set keep to 0 for all trials
        keepStruct(trial).keep = 0;
        %Trial Type
        trialName = Data(trial).trialName;
        for i = 1:numTrialTypes
           if strcmpi(trialInclStates(i).trialName,trialName)
               trialType = i;
               keepStruct(trial).trialType = trialType;
           end
        end
        %Length
        keepStruct(trial).length = 0;
        inclStates = trialInclStates(trialType).inclStates;
        numInclStates = size(inclStates,2);
        stateTransitions = double(Data(trial).stateData.stateTransitions);
        stateNames = Data(trial).stateData.stateNames;
        for state = 1:numInclStates
            stateInd = max(find([cellfun(@(x) strcmpi(x,inclStates{1,state}),stateNames)]==1));
            occurrence = trialInclStates(trialType).inclOccurrence{1,state};
            if ismember(stateInd,stateTransitions(1,:))
                %Get State Start and End Times 
                if strcmpi(occurrence,'first')
                    stateStartTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==stateInd))));
                    stateEndTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==stateInd))+1));
                elseif strcmpi(occurrence,'last')
                    stateStartTime = double(stateTransitions(2,max(find(stateTransitions(1,:)==stateInd))));
                    stateEndTime = double(stateTransitions(2,max(find(stateTransitions(1,:)==stateInd))+1));
                end
            end
            keepStruct(trial).length = (stateEndTime-stateStartTime) + keepStruct(trial).length;
        end
    end
    
    %Determine trials to keep for each trial type
    for trialType = 1:numTrialTypes
        tempData = keepStruct([keepStruct.trialType]==trialType);
        lengths = [tempData.length];
        meanLength = mean(lengths);
        stdLength = std(lengths);
        minKeepLength = meanLength - 2*stdLength;
        if minKeepLength < 0
            minKeepLength = 0;
        end
        maxKeepLength = meanLength + 2*stdLength;
        %Plot distribution and keep lengths
        figure
            histogram(lengths)
            hold on
            ca = gca;
            maxY = ca.YLim(2);
            line([minKeepLength,minKeepLength],[0 maxY],'Color','red')
            line([maxKeepLength,maxKeepLength],[0 maxY],'Color','red')
            xlabel('Length (ms)')
            ylabel('Count');
            title(['trialType ',num2str(trialType),' Lengths'])
        for trial = 1:numTrials
           curTrialType = keepStruct(trial).trialType;
           length = keepStruct(trial).length;
           if curTrialType == trialType & length > minKeepLength & length < maxKeepLength
               keepStruct(trial).keep = 1;
           end
        end
    end

    %Remove trials
    Data = Data([keepStruct.keep]==1);

end