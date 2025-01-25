function [sortMean,sortStd,sortRange] = getZScoreParams(Data,binWidth,kernelStdDev)

%This function smooths and downsamples neural data to get estimates of
%firing rates for each sort. It then gets means and stdDevs for z-scoring
%and removal of low-firing channels.

%% Set up trialInclStates
    trialInclStates = struct('trialName','','inclStates',[]);
    
%% Get all smoothed neural data
    allNeuralData = [];
    numTrials = size(Data,2);
    for trial = 1:numTrials
        trial
        %Get trial name
        trialName = Data(trial).trialName;
        trialInclStates(1).trialName = {trialName};
        %Match trial name with appropriate inclStates
        switch trialName
            %Earl
            case 'GridTask_CO_Across_BC_ForceBar'
                trialInclStates(1).inclStates = {{'state','Center Target','first',0},{'state','Success with Reward','first',0}};
                %Tried changing to included more states on 10/31/23 but
                %there wasn't enought time in trials to z-score if I did
                %this. These are commented windows below.
                %trialInclStates(1).inclStates = {{'state','Touch Bar Hold','first',0},{'state','Success with Reward','first',0}};
            case 'GridReaching'
                trialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Success with Reward','first',0}};
                %trialInclStates(1).inclStates = {{'state','Center Acquire','first',0},{'state','Success with Reward','first',0}};
            case 'GridTask_BC_ForceBar' %From multiple tasks paradigm
                trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',0},{'state','Success with Reward','first',0}};
                %trialInclStates(1).inclStates = {{'state','Touch Bar Hold','first',0},{'state','Success with Reward','first',0}};
            case 'HC_CenterOut_ForceBar_20200314' %From multiple tasks paradigm
                trialInclStates(1).inclStates = {{'state','Touch Bar Hold','first',0},{'state','Success with Reward','first',0}};              
            case 'IsometricForce_1D' %From multiple tasks paradigm
                trialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Success with Reward','first',0}};
                %trialInclStates(1).inclStates = {{'state','Center','first',0},{'state','Success with Reward','first',0}};
            case 'BCI Center Out' %From shoulder/elbow BCI
                trialInclStates(1).inclStates = {{'state','Center Target','first',0},{'state','Success with Reward','first',0}};
                %trialInclStates(1).inclStates = {{'state','Touch Bar Hold','first',0},{'state','Success with Reward','first',0}};
            case 'DelayedCenterOut20210828' %From shoulder/elbow BCI + DCO session
                trialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Success with Reward','first',0}};
                %trialInclStates(1).inclStates = {{'state','Center Acquire','first',0},{'state','Success with Reward','first',0}};
            case 'CenterOut_ForceBar_BC' %Earl CO BCI Used for drift control (20190830)
                trialInclStates(1).inclStates = {{'state','BC Freeze','first',0},{'state','Success with Reward','first',0}};   
                %trialInclStates(1).inclStates = {{'state','Touch Bar Hold','first',0},{'state','Success with Reward','first',0}};  
            case 'CentetOut_BC_TouchBar' %Earl CO BCI Used for for alternate decoder exps
                trialInclStates(1).inclStates = {{'state','Center Target','first',0},{'state','Success with Reward','first',0}}; 
                %trialInclStates(1).inclStates = {{'state','Touch Bar Hold','first',0},{'state','Success with Reward','first',0}}; 
            case 'CenterOutCenter_BC_TouchBar' %Earl COC BCI Used for for alternate decoder exps
                trialInclStates(1).inclStates = {{'state','Center Target','first',0},{'state','Success with Reward','first',0}};     
                %trialInclStates(1).inclStates = {{'state','Touch Bar Hold','first',0},{'state','Success with Reward','first',0}}; 
            %Nigel
            case 'Nigel Dissociation'
                trialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Success with Reward','first',0}};
                %trialInclStates(1).inclStates = {{'state','Trial Start','first',0},{'state','Success with Reward','first',0}};
            case 'Nigel Posture BC Center Out'
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Success with Reward','first',0}};
                %trialInclStates(1).inclStates = {{'state','Start','first',0},{'state','Success with Reward','first',0}};
            %Rocky
            case 'Rocky Dissociation'
                trialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Success with Reward','first',0}};
                %trialInclStates(1).inclStates = {{'state','Trial Start','first',0},{'state','Success with Reward','first',0}};
            case 'Rocky Posture BC Center Out'
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
                %trialInclStates(1).inclStates = {{'state','Hold A','first',0},{'state','Hold','first',0}};
            %No case found
            otherwise
                error('No case for this trial name in getZScoreParams')
                
        end
        %Get associated neural data
        [statesTraj,~] = getStatesTraj20220419(Data(trial),trialInclStates,'smoothFR',binWidth,kernelStdDev);
        allNeuralData = vertcat(allNeuralData,statesTraj);
    end
    
%% Get all params
    sortMean = mean(allNeuralData,1);
    sortStd = std(allNeuralData,1);
    sortRange = max(allNeuralData) - min(allNeuralData);
    
end