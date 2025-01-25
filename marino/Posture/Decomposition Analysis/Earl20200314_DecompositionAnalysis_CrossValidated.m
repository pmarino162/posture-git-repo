clear; clc; clf; close all;

%% Load, Preprocess, and Label Data    
    %3/14/2020
    date = '20200314';
    exclCh =  [3 12 16 20 31 73 85 92 94];
    [Data,I30GPFAParams,I30DecoderParams,...
     N00GPFAParams,N00DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200314;
    numPostures = 5;
    numTargets = 8;

%% Prediction 
    numFolds = 5;
    for taskID = 1:3
        taskData = Data([Data.taskData.taskID]==taskID);
        numObs = size(taskData,1);
        partition = cvpartition(numObs,'KFold',numFolds);
        for fold = 1:numFolds
           %Get training and test sets
           trainInds = training(partition,fold); testInds = test(partition,fold);
           trainData = taskData(trainInds); testData = taskData(testInds);
           %Get condition averages
           switch taskID
               case 1
                    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                    trialInclStates(1).inclStates = {'Step 1'}; 
                    trialInclStates(1).inclOccurrence = {'last'};
               case 2
                    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
                    trialInclStates(2).inclStates = {'Reach'};
                    trialInclStates(2).inclOccurrence = {'first'};
               case 3
                    trialInclStates(3).trialName = {'IsometricForce_1D'};
                    trialInclStates(3).inclStates = {'Target'};
                    trialInclStates(3).inclOccurrence = {'last'};
           end
           condFields = {{'posture','postureData','postureID'},{'target','targetData','targetID'}};
           trajFields = {'allChannelSmoothedFR'};
           trialInclStates(1).addTimeToBeginning = {-100};
           trialInclStates(1).addTimeToEnd = {0};
           trainTrajStruct = getTrajStruct(trainData,condFields,trajFields,trialInclStates);
           testTrajStruct = getTrajStruct(testData,condFields,trajFields,trialInclStates);
           %Get postural offset
           
           %Get CI component
           
           %Get target component
           
           %Create prediction for each condition
           
           %Assess model performance
           
           %Save results
        end
    end


%% Plotting 
    
