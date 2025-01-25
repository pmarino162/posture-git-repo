function [Data] = getKinData20220419(Data)

%For every trial in data struct, get kinematic features of interest 
%(reaction time, peak speed time, etc). Add to
%Data struct in kinData field.


numTrials = size(Data,2);
binWidth = 1; %Necessary input to getStatesTraj20220419, but not used
kernelStdDev = NaN; %Necessary input to getStatesTraj20220419, but not used
rmTrials = [];

for trial = 1:numTrials
    kinData = [];
    trialName = Data(trial).trialName;
    
    %Get state data and set up trialInclStates
    stateNames = Data(trial).stateData.stateNames;
    stateTransitions = Data(trial).stateData.stateTransitions;
    trialInclStates = struct('trialName','','inclStates',[]);
    trialInclStates(1).trialName = {trialName};
    
    switch trialName
        %BCI for all monkeys, incl. Earl multiple-task 
        case {'GridTask_CO_Across_BC_ForceBar','GridTask_BC_ForceBar','Nigel Posture BC Center Out','Rocky Posture BC Center Out',...
                'BCI Center Out','CentetOut_BC_TouchBar','CenterOutCenter_BC_TouchBar','CenterOut_ForceBar_BC'}
            %Get trialType specific info
            if strcmp(trialName,'GridTask_CO_Across_BC_ForceBar')
               movePeriodStartID = max(find([cellfun(@(x) strcmpi(x,'Step 1'),stateNames)]==1));
               movePeriodEndID = max(find([cellfun(@(x) strcmpi(x,'Step 2'),stateNames)]==1));
            end
            if strcmp(trialName,'GridTask_BC_ForceBar') %Earl multiple tasks BCI
               movePeriodStartID = max(find([cellfun(@(x) strcmpi(x,'Step 1 Freeze'),stateNames)]==1));
               movePeriodEndID = max(find([cellfun(@(x) strcmpi(x,'Step 2'),stateNames)]==1));
            end
            if strcmp(trialName,'Nigel Posture BC Center Out')
               movePeriodStartID = max(find([cellfun(@(x) strcmpi(x,'Cursor Release'),stateNames)]==1));
               movePeriodEndID = max(find([cellfun(@(x) strcmpi(x,'Target Hold'),stateNames)]==1));
            end
            if strcmp(trialName,'Rocky Posture BC Center Out')
               movePeriodStartID = max(find([cellfun(@(x) strcmpi(x,'React'),stateNames)]==1));
               movePeriodEndID = max(find([cellfun(@(x) strcmpi(x,'Hold'),stateNames)]==1));
            end
            if strcmp(trialName,'BCI Center Out')
               movePeriodStartID = max(find([cellfun(@(x) strcmpi(x,'Step 1'),stateNames)]==1));
               movePeriodEndID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
            end
            if strcmp(trialName,'CentetOut_BC_TouchBar')
               movePeriodStartID = max(find([cellfun(@(x) strcmpi(x,'Step 1'),stateNames)]==1));
               movePeriodEndID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
            end
            if strcmp(trialName,'CenterOutCenter_BC_TouchBar')
               movePeriodStartID = max(find([cellfun(@(x) strcmpi(x,'Step 1'),stateNames)]==1));
               movePeriodEndID = max(find([cellfun(@(x) strcmpi(x,'Step 2'),stateNames)]==1));
            end
            if strcmp(trialName,'CenterOut_ForceBar_BC')
               movePeriodStartID = max(find([cellfun(@(x) strcmpi(x,'BC Freeze'),stateNames)]==1));
               movePeriodEndID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
            end
            movePeriodStartTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==movePeriodStartID))));
            movePeriodEndTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==movePeriodEndID))));
            %Store results in kinData
            kinData.moveTime = movePeriodEndTime - movePeriodStartTime;
        
        %Earl Isometric force
        case {'IsometricForce_1D'}
            %Get event times 
            centerHoldID = max(find([cellfun(@(x) strcmpi(x,'Center Hold'),stateNames)]==1));
            targetAcqID = max(find([cellfun(@(x) strcmpi(x,'Target'),stateNames)]==1));
            targetHoldID = max(find([cellfun(@(x) strcmpi(x,'Target Hold'),stateNames)]==1));
            successID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));         
            
            %Get event times 
            centerEnterTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==centerHoldID))));
            goTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetAcqID))));
            targetEnterTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetHoldID))));
            successTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==successID))));
            
            %Get kinematic info from reach, list of ill-behaved trials to remove
            [moveOnsetTime,rmTrialFlag] = getIsoKin(Data(trial),binWidth,kernelStdDev);   
            if rmTrialFlag
                rmTrials = [rmTrials,trial];
            end
            
            %Store results in kinData
            kinData.centerHoldTime = goTime - centerEnterTime;
            kinData.targetHoldTime = successTime - targetEnterTime;
            kinData.moveOnsetTime = moveOnsetTime;
            kinData.rxnTime = moveOnsetTime - goTime;
            kinData.reachTime = targetEnterTime - moveOnsetTime;
            kinData.moveTime = targetEnterTime - goTime;
        
        %Nigel Dissociation; Earl 7-posture DCO; Earl 3-target reaching in
        %multiple tasks paradigm
        case {'Nigel Dissociation','Rocky Dissociation','GridReaching','HC_CenterOut_ForceBar_20200314'}
            %Get trialType specific info
            if strcmp(trialName,'Nigel Dissociation') | strcmp(trialName,'Rocky Dissociation')
                centerHoldID = max(find([cellfun(@(x) strcmpi(x,'Center Hold'),stateNames)]==1));
                targetAcqID = max(find([cellfun(@(x) strcmpi(x,'Reach'),stateNames)]==1));
                targetHoldID = max(find([cellfun(@(x) strcmpi(x,'Target Hold'),stateNames)]==1));
                successID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
                trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Target Hold','first',0}};
            elseif strcmp(trialName,'GridReaching')
                centerHoldID = max(find([cellfun(@(x) strcmpi(x,'Center Hold'),stateNames)]==1));
                delayID = max(find([cellfun(@(x) strcmpi(x,'Delay'),stateNames)]==1));
                targetAcqID = max(find([cellfun(@(x) strcmpi(x,'Target Acquire'),stateNames)]==1));
                targetHoldID = max(find([cellfun(@(x) strcmpi(x,'Target Hold'),stateNames)]==1));
                successID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
                trialInclStates(1).inclStates = {{'state','Target Acquire','first',0},{'state','Target Hold','first',0}};
            elseif strcmp(trialName,'HC_CenterOut_ForceBar_20200314')
                centerHoldID = max(find([cellfun(@(x) strcmpi(x,'Touch Bar Hold'),stateNames)]==1));
                targetAcqID = max(find([cellfun(@(x) strcmpi(x,'Reach'),stateNames)]==1));
                targetHoldID = max(find([cellfun(@(x) strcmpi(x,'Hold'),stateNames)]==1));
                successID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
                trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Hold','first',0}};
            end
            %Get event times
            if strcmp(trialName,'GridReaching')
               delayTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==delayID))));
            end
            centerEnterTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==centerHoldID))));
            goTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetAcqID))));
            targetEnterTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetHoldID))));
            successTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==successID))));
            %Get kinematic info from reach, list of ill-behaved trials to remove
            [moveOnsetTime,peakSpeedTime,rmTrialFlag] = getReachKin(Data(trial),trialInclStates,binWidth,kernelStdDev);   
            if rmTrialFlag
                rmTrials = [rmTrials,trial];
            end
            %Store results in kinData
            if strcmp(trialName,'GridReaching')
               kinData.delayLength = goTime-delayTime;
               kinData.centerHoldTime = delayTime - centerEnterTime;
            elseif strcmp(trialName,'Nigel Dissociation') || strcmp(trialName,'Rocky Dissociation') || strcmp(trialName,'HC_CenterOut_ForceBar_20200314')
               kinData.centerHoldTime = goTime - centerEnterTime;
            end
            kinData.targetHoldTime = successTime - targetEnterTime;
            kinData.moveOnsetTime = moveOnsetTime;
            kinData.rxnTime = moveOnsetTime - goTime;
            kinData.reachTime = targetEnterTime - moveOnsetTime;
            kinData.peakSpeedTime = peakSpeedTime;
    end
    Data(trial).kinData = kinData;
end

%% Remove ill-behaved trials
    Data(rmTrials) = [];
    
end