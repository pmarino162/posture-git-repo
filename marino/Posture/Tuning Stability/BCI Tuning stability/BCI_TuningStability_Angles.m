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
    [Data] = subselectForTrajDist(Data,task);
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

%% Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
     
%% Create struct with tuning axes and varExpl
    %Target tuning across postures
    targetAxesStruct = struct('posture',[],'targetMeans',[],'axes',[],'explained',[]);
    structInd = 1;
    for posture = postureList
         tempData = trajStruct([trajStruct.posture]==posture);
         allAvgs = [];
         for i = 1:size(tempData,2)
             allAvgs = vertcat(allAvgs,tempData(i).avgSmoothFR.condAvg);
         end
         [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs); 
         targetAxesStruct(structInd).posture = posture;
         targetAxesStruct(structInd).targetMeans = allAvgs;
         targetAxesStruct(structInd).axes = coeff;
         targetAxesStruct(structInd).explained = explained;
         structInd = structInd + 1;
    end

    %Posture tuning across targets
    postureAxesStruct = struct('target',[],'postureMeans',[],'axes',[],'explained',[]);
    structInd = 1;
    for target = targetList
         tempData = trajStruct([trajStruct.target]==target);
         allAvgs = [];
         for i = 1:size(tempData,2)
             allAvgs = vertcat(allAvgs,tempData(i).avgSmoothFR.condAvg);
         end
         [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs); 
         postureAxesStruct(structInd).target = target;
         postureAxesStruct(structInd).postureMeans = allAvgs;
         postureAxesStruct(structInd).axes = coeff;
         postureAxesStruct(structInd).explained = explained;
         structInd = structInd + 1;
    end
    
    
    
    %% Create obsStruct for classification
    obsStruct = struct('target',[],'posture',[],'postureLabelVec',[],'targetLabelVec',[],'numObs',[],'allObs',[]);
    structInd = 1; 
    for target = targetList
        for posture = postureList
            obsStruct(structInd).target = target;
            obsStruct(structInd).posture = posture;
            allObs = [];
            tempTrajStruct = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture);
            for i = 1:size(tempTrajStruct,2)
                for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
                    trialAvg = tempTrajStruct(i).allSmoothFR(j).trialAvg;
                    allObs = vertcat(allObs,trialAvg);
                end
            end
            obsStruct(structInd).allObs = allObs;
            numObs = size(obsStruct(structInd).allObs,1);
            obsStruct(structInd).numObs = numObs;
            obsStruct(structInd).postureLabelVec = ones(numObs,1)*posture;
            obsStruct(structInd).targetLabelVec = ones(numObs,1)*target;
            structInd = structInd + 1;
        end
    end

%% Train posture classifier on individual target; classify across targets
    acrossTargetResult = zeros(numTargets,numTargets);
    for trainTarget = targetList
       %Train model on current target
       observations = vertcat(obsStruct([obsStruct.target]==trainTarget).allObs);
       labels = vertcat(obsStruct([obsStruct.target]==trainTarget).postureLabelVec);
       [LDAModel,numClasses] = LDATrainModel(observations,labels);
       %Test across all other targets
       for testTarget = targetList
          observations = []; labels = [];
          minNumObs = min([obsStruct([obsStruct.target]==testTarget).numObs]);
          for posture = postureList
            numObs = obsStruct([obsStruct.target]==testTarget & [obsStruct.posture]==posture).numObs;
            observations = vertcat(observations,obsStruct([obsStruct.target]==testTarget & [obsStruct.posture]==posture).allObs(randsample(numObs,minNumObs),:));
            labels = vertcat(labels,obsStruct([obsStruct.target]==testTarget & [obsStruct.posture]==posture).postureLabelVec(1:minNumObs));
          end
          [predictedLabels,pctCorrect] = LDAClassify(observations,labels,numClasses,LDAModel);
          acrossTargetResult(trainTarget,testTarget) = pctCorrect;
       end
    end


%% Train target classifier on individual postures; classify across postures
    acrossPostureResult = zeros(numPostures,numPostures);
    for trainPosture = postureList
       %Train model on current posture
       observations = vertcat(obsStruct([obsStruct.posture]==trainPosture).allObs);
       labels = vertcat(obsStruct([obsStruct.posture]==trainPosture).targetLabelVec);
       [LDAModel,numClasses] = LDATrainModel(observations,labels);
       %Test across all other postures
       for testPosture = postureList
          observations = []; labels = [];
          minNumObs = min([obsStruct([obsStruct.posture]==testPosture).numObs]);
          for target = targetList
            numObs = obsStruct([obsStruct.posture]==testPosture & [obsStruct.target]==target).numObs;
            observations = vertcat(observations,obsStruct([obsStruct.posture]==testPosture & [obsStruct.target]==target).allObs(randsample(numObs,minNumObs),:));
            labels = vertcat(labels,obsStruct([obsStruct.posture]==testPosture & [obsStruct.target]==target).postureLabelVec(1:minNumObs));
          end
          [predictedLabels,pctCorrect] = LDAClassify(observations,labels,numClasses,LDAModel);
          acrossPostureResult(trainPosture,testPosture) = pctCorrect;
       end
    end

    
%% Plot across-target posture classification performance    
    figure
    clims = [20,100];
    colormap winter
    imagesc(acrossTargetResult,clims)
    xticks([1:numTargets])
    xlabel('Testing Target')
    yticks([1:numTargets])
    ylabel('Training Target')
    c = colorbar;
    c.Ticks = [20 50 75 100];
%     c.TickLabels = {'Chance (33.33)','50','75','100'};
    c.Label.String = 'Classification Accuracy (%)';
    
%% Plot across-posture target classification performance    
    figure
    clims = [12.5,100];
    colormap winter
    imagesc(acrossPostureResult,clims)
    xticks([1:numPostures])
    xlabel('Testing Posture')
    yticks([1:numPostures])
    ylabel('Training Posture')
    c = colorbar;
    c.Ticks = [12.5 50 75 100];
%     c.TickLabels = {'Chance (33.33)','50','75','100'};
    c.Label.String = 'Classification Accuracy (%)';