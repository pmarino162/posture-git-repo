clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;

%% Load Data, subselect
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;
    task = 'BCI';
    [Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
%     [Data] = subselectForTrajDist(Data,task);
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
    %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-150},{'state','Step 2','first',0}};
    trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%% For each timecourse in trajstruct, get single point
    for i = 1:size(trajStruct,2)
        %Individual trials
        for j = 1:size(trajStruct(i).allSmoothFR,2)
           trajStruct(i).allSmoothFR(j).trialAvg = mean(trajStruct(i).allSmoothFR(j).traj); 
        end
        %Condition average
        trajStruct(i).avgSmoothFR.condAvg = mean(trajStruct(i).avgSmoothFR.traj);
    end
    numDims = size(trajStruct(1).avgSmoothFR.condAvg,2);
    
%% Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
     
%% Create obsStruct for classification
    obsStruct = struct('target',[],'posture',[],'postureLabelVec',[],'targetLabelVec',[],'numObs',[],'allObs',[]);
    structInd = 1; 
    for target = targetList
        for posture = postureList
            obsStruct(structInd).target = target;
            obsStruct(structInd).posture = posture;
            allObs = [];
            tempTrajStruct = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture);
            for i = 1:size(tempTrajStruct.allSmoothFR,2)
                trialAvg = tempTrajStruct.allSmoothFR(i).trialAvg;
                allObs = vertcat(allObs,trialAvg);
            end
            obsStruct(structInd).allObs = allObs;
            numObs = size(obsStruct(structInd).allObs,1);
            obsStruct(structInd).numObs = numObs;
            obsStruct(structInd).postureLabelVec = ones(numObs,1)*posture;
            obsStruct(structInd).targetLabelVec = ones(numObs,1)*target;
            structInd = structInd + 1;
        end
    end

%% Get target-specific train and test sets for posture classification
    targetObsStruct = struct('target',[],'obs',[],'labels',[],'partition',[]);
    numFolds = 5;
    structInd = 1;
    for target = targetList
        minNumObs = min([obsStruct([obsStruct.target]==target).numObs]);
        obs = NaN(minNumObs*numPostures,numDims);
        labels = NaN(minNumObs*numPostures,1);
        i = 1;
        for posture = postureList
            numObs = obsStruct([obsStruct.target]==target & [obsStruct.posture]==posture).numObs;
            obs(i:i+minNumObs-1,:) = vertcat(obsStruct([obsStruct.target]==target & [obsStruct.posture]==posture).allObs(randsample(numObs,minNumObs),:));
            labels(i:i+minNumObs-1,1) = vertcat(obsStruct([obsStruct.target]==target & [obsStruct.posture]==posture).postureLabelVec(1:minNumObs));
            i = i + minNumObs;
        end
        c = cvpartition(minNumObs*numPostures,'KFold',numFolds);
        targetObsStruct(structInd).target = target;
        targetObsStruct(structInd).obs = obs;
        targetObsStruct(structInd).labels = labels;
        targetObsStruct(structInd).partition = c;
        structInd = structInd + 1;
    end

%% Get posture-specific train and test sets for target classification
    postureObsStruct = struct('posture',[],'obs',[],'labels',[],'partition',[]);
    numFolds = 5;
    structInd = 1;
    for posture = postureList
        minNumObs = min([obsStruct([obsStruct.posture]==posture).numObs]);
        obs = NaN(minNumObs*numTargets,numDims);
        labels = NaN(minNumObs*numTargets,1);
        i = 1;
        for target = targetList
            numObs = obsStruct([obsStruct.target]==target & [obsStruct.posture]==posture).numObs;
            obs(i:i+minNumObs-1,:) = vertcat(obsStruct([obsStruct.target]==target & [obsStruct.posture]==posture).allObs(randsample(numObs,minNumObs),:));
            labels(i:i+minNumObs-1,1) = vertcat(obsStruct([obsStruct.target]==target & [obsStruct.posture]==posture).targetLabelVec(1:minNumObs));
            i = i + minNumObs;
        end
        c = cvpartition(minNumObs*numTargets,'KFold',numFolds);
        postureObsStruct(structInd).posture = posture;
        postureObsStruct(structInd).obs = obs;
        postureObsStruct(structInd).labels = labels;
        postureObsStruct(structInd).partition = c;
        structInd = structInd + 1;
    end
    
%% Train posture classifier on individual target; classify across targets
    acrossTargetResult = zeros(numTargets,numTargets,numFolds);
    for trainTarget = targetList
       trainObs = targetObsStruct([targetObsStruct.target]==trainTarget).obs;
       trainLabels = targetObsStruct([targetObsStruct.target]==trainTarget).labels;
       trainPartition = targetObsStruct([targetObsStruct.target]==trainTarget).partition;
       for testTarget = targetList
            testObs = targetObsStruct([targetObsStruct.target]==testTarget).obs;
            testLabels = targetObsStruct([targetObsStruct.target]==testTarget).labels;
            testPartition = targetObsStruct([targetObsStruct.target]==testTarget).partition;
            for fold = 1:numFolds
               [NBModel,numClasses] = NBTrainModel(trainObs(trainPartition.training(fold),:),trainLabels(trainPartition.training(fold),:)); 
               [predictedLabels,pctCorrect] = NBClassify(testObs(testPartition.test(fold),:),testLabels(testPartition.test(fold),:),numClasses,NBModel);
               acrossTargetResult(trainTarget,testTarget,fold) = pctCorrect;
            end   
       end
    end    
    acrossTargetResult = mean(acrossTargetResult,3);

%% Train target classifier on individual posture; classify across postures
    acrossPostureResult = zeros(numPostures,numPostures,numFolds);
    for trainPosture = postureList
       trainObs = postureObsStruct([postureObsStruct.posture]==trainPosture).obs;
       trainLabels = postureObsStruct([postureObsStruct.posture]==trainPosture).labels;
       trainPartition = postureObsStruct([postureObsStruct.posture]==trainPosture).partition;
       for testPosture = postureList
            testObs = postureObsStruct([postureObsStruct.posture]==testPosture).obs;
            testLabels = postureObsStruct([postureObsStruct.posture]==testPosture).labels;
            testPartition = postureObsStruct([postureObsStruct.posture]==testPosture).partition;
            for fold = 1:numFolds
               [NBModel,numClasses] = NBTrainModel(trainObs(trainPartition.training(fold),:),trainLabels(trainPartition.training(fold),:)); 
               [predictedLabels,pctCorrect] = NBClassify(testObs(testPartition.test(fold),:),testLabels(testPartition.test(fold),:),numClasses,NBModel);
               acrossPostureResult(trainTarget,testTarget,fold) = pctCorrect;
            end   
       end
    end    
    acrossPostureResult = mean(acrossPostureResult,3);
    
%% Plot across-target posture classification performance    
    figure
    clims = [20,max(max(acrossTargetResult))];
    colormap hot
    imagesc(acrossTargetResult,clims)
    xticks([1:numTargets])
    xlabel('Testing Target')
    yticks([1:numTargets])
    ylabel('Training Target')
    c = colorbar;
    c.Label.String = 'Classification Accuracy (%)';
    ax = gca;
    ax.FontSize = 16;
    
%% Plot across-posture target classification performance    
    figure
    clims = [12.5,max(max(acrossPostureResult))];
    colormap hot
    imagesc(acrossPostureResult,clims)
    xticks([1:numPostures])
    xlabel('Testing Posture')
    yticks([1:numPostures])
    ylabel('Training Posture')
    c = colorbar;
    c.Label.String = 'Classification Accuracy (%)';