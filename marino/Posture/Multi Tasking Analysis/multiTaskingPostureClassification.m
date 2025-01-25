clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;

%% Load Data, subselect
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;
    
    [Data,I30GPFAParams,I30DecoderParams,...
     N00GPFAParams,N00DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200314_20211210;

%% Get trajStruct
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    
    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
        trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',0},{'state','Step 2','first',0}};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(2).inclStates = {{'state','Reach','first',0},{'state','Hold','first',0}};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
    
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


%% Get posture list
    postureList = unique([trajStruct.posture]);

%% Do LDA by posture on all trial averages 
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = postureList
       obsStruct(structInd).label = posture;
       allObs = [];
       tempTrajStruct = trajStruct([trajStruct.posture]==posture);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               trialAvg = tempTrajStruct(i).allSmoothFR(j).trialAvg;
               allObs = vertcat(allObs,trialAvg);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    [postureLDA] = doLDA(obsStruct);
    
    
%% Plot in LDA Space   
    figure
    hold on
    for i = 1:size(trajStruct,2)
        posture = trajStruct(i).posture;
        task = trajStruct(i).task;
        allTrialAvg = vertcat(trajStruct(i).allSmoothFR.trialAvg);
        ldaProj = allTrialAvg*postureLDA;
        if task == 1
            plot(ldaProj(:,1),ldaProj(:,2),'s','MarkerSize',7,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
        elseif task == 2
            plot(ldaProj(:,1),ldaProj(:,2),'>','MarkerSize',7,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
        elseif task == 3 
            plot(ldaProj(:,1),ldaProj(:,2),'p','MarkerSize',7,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
        end
    end 
    xlabel('Posture LDA 1')
    ylabel('Posture LDA 2')
    
%% Create obsStruct for across task classification
%NEED TO BALANCE ACROSS TARGETS
    taskObsStruct = struct('task',[],'posture',[],'postureLabelVec',[],'numObs',[],'allObs',[]);
    structInd = 1; 
    for task = 1:3
        for posture = postureList
            taskObsStruct(structInd).task = task;
            taskObsStruct(structInd).posture = posture;
            allObs = [];
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            for i = 1:size(tempTrajStruct,2)
                for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
                    trialAvg = tempTrajStruct(i).allSmoothFR(j).trialAvg;
                    allObs = vertcat(allObs,trialAvg);
                end
            end
            taskObsStruct(structInd).allObs = allObs;
            numObs = size(taskObsStruct(structInd).allObs,1);
            taskObsStruct(structInd).numObs = numObs;
            taskObsStruct(structInd).postureLabelVec = ones(numObs,1)*posture;
            structInd = structInd + 1;
        end
    end


%% Train on individual tasks, classify across tasks
    resultStruct = struct('trainTask',[],'testTask',[],'pctCorrect',[]);
    resultMat = zeros(3,3);
    structInd = 1;
    for trainTask = 1:3
       %Train model on current task
       observations = vertcat(taskObsStruct([taskObsStruct.task]==trainTask).allObs);
       labels = vertcat(taskObsStruct([taskObsStruct.task]==trainTask).postureLabelVec);
       [LDAModel,numClasses] = LDATrainModel(observations,labels);
       %Test on each task individually
       for testTask = 1:3
          observations = vertcat(taskObsStruct([taskObsStruct.task]==testTask).allObs);
          labels = vertcat(taskObsStruct([taskObsStruct.task]==testTask).postureLabelVec);
          [predictedLabels,pctCorrect] = LDAClassify(observations,labels,numClasses,LDAModel);
          resultMat(trainTask,testTask) = pctCorrect;
          resultStruct(structInd).trainTask = trainTask;
          resultStruct(structInd).testTask = testTask;
          resultStruct(structInd).pctCorrect = pctCorrect;
          structInd = structInd + 1; 
       end
        
    end
    
%% Plot classification performance    
    figure
    clims = [33.33,100];
    colormap winter
    imagesc(resultMat,clims)
    xticks([1:3])
    xticklabels({'BCI','Reaching','Isometric Force'})
    xlabel('Testing Set')
    yticks([1:3])
    yticklabels({'BCI','Reaching','Isometric Force'})
    ylabel('Training Set')
    c = colorbar;
    
    c.Ticks = [33.33 50 75 100];
%     c.TickLabels = {'Chance (33.33)','50','75','100'};
    c.Label.String = 'Classification Accuracy (%)';
    
%% Local function for performing LDA
     function [LDAproj] = doLDA(obsStruct)
        %Preallocate
        minNumObs = min([obsStruct.numObs]);
        numClasses = size(obsStruct,2);
        numDims = size(obsStruct(1).allObs,2);
        obs = NaN(minNumObs*numClasses,numDims); labels =  NaN(minNumObs*numClasses,1); 
        %Fill
        k = 1;
        for i = 1:size(obsStruct,2)
            totalClassObs = size(obsStruct(i).allObs,1);
            obs(k:k+minNumObs-1,:) = obsStruct(i).allObs(randsample(totalClassObs,minNumObs),:);
            labels(k:k+minNumObs-1,1) = ones(minNumObs,1).*obsStruct(i).label;
            k = k+minNumObs;
        end
        LDAproj = fisherLDA(obs, labels);
        [LDAproj,~] = qr(LDAproj);
        LDAproj = LDAproj(:,1:numClasses-1);
     end  