clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Error Reduction Figures';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat');
    pcmap = orli;
    
%% Parameters
   numManPCs = 10;     %Num PCs to project data into before analysis
   numFolds = 5;
   
%% Setup resultStruct
    resultStruct = struct('dataset',[],'result',[]);
    resultStructInd = 1;
    
%% Load and preprocess data
    for datasetList = {'E20200314'}%{'E20200311','E20200312','E20200313','E20200314'}
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25; kernelStdDev = 25;
        trialInclStates = struct('trialName','','inclStates',[]);
        trajFields = {'zSmoothFR','markerVel','marker'};
        condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};

        trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
            trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',0},{'state','Step 2','first',0}};
        trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
            trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-300},{'kin','moveOnsetTime','first',-100}};
        trialInclStates(3).trialName = {'IsometricForce_1D'};
            trialInclStates(3).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

        taskIDs = struct('ID',[],'task','');
        taskIDs(1).ID = 1; taskIDs(1).task = 'BC';
        taskIDs(2).ID = 2; taskIDs(2).task = 'HC';
        taskIDs(3).ID = 3; taskIDs(3).task = 'Iso';
        taskIDs(4).ID = 4; taskIDs(4).task = 'All';

        %For each timecourse, get one point
         for i = 1:size(trajStruct,2)
            %Individual trials
            for j = 1:size(trajStruct(i).allSmoothFR,2)
               trajStruct(i).allSmoothFR(j).obs = mean(trajStruct(i).allSmoothFR(j).traj); 
            end
            %Condition average
            trajStruct(i).avgSmoothFR.condAvg = mean(trajStruct(i).avgSmoothFR.traj);
         end

        %Reduce dimensionality of observations
        allObs = [];
        for i = 1:size(trajStruct,2)
            for j = 1:size(trajStruct(i).allSmoothFR,2)
                obs = trajStruct(i).allSmoothFR(j).obs;
                allObs = vertcat(allObs,obs);
            end
        end
        [coeff,score,latent,tsquared,explained,mu] = pca(allObs); 
        for i = 1:size(trajStruct,2)
            for j = 1:size(trajStruct(i).allSmoothFR,2)
                trajStruct(i).allSmoothFR(j).obs = trajStruct(i).allSmoothFR(j).obs*coeff(:,1:numManPCs);
            end
        end

        %Get posture list
        postureList = unique([trajStruct.posture]);
    
        %% Create taskObsStruct, balancing observations by task when task=4 (all tasks)
        %Get mininum number of trials for task*posture condition
        minNumCondObs = zeros(1,3);
        numCond = zeros(1,3);
        numPossibleObsToUse = zeros(1,3);
        for task = 1:3
            tempMinNumCondObs = 999999999;
            tempTrajStruct = trajStruct([trajStruct.task]==task);
            for i = 1:size(tempTrajStruct,2)
                numCondObs = size(vertcat(tempTrajStruct(i).allSmoothFR.obs),1);
                if numCondObs < tempMinNumCondObs
                    tempMinNumCondObs = numCondObs;
                end
            end
            minNumCondObs(task) = tempMinNumCondObs;
            numCond(task) = size(tempTrajStruct,2);
            numPossibleObsToUse(task) = minNumCondObs(task)*numCond(task);
        end    
        numObsToUse = min(numPossibleObsToUse);

        %Build taskObsStruct
        taskObsStruct = struct('task',[],'obs',[],'labels',[],'c',[]);
        structInd = 1;
        for task = 1:4
           obs = []; labels = [];
           if task == 4
                tempTrajStruct = trajStruct;
                numCond = size(tempTrajStruct,2);
                numObsPerCond = floor(numObsToUse/numCond);
                for i = 1:size(tempTrajStruct,2)
                    condObs = vertcat(tempTrajStruct(i).allSmoothFR.obs);
                    obs = vertcat(obs,datasample(condObs,numObsPerCond,'Replace',false));
                    posture = tempTrajStruct(i).posture;
                    labels = vertcat(labels,ones(numObsPerCond,1)*posture);
                end      
           else
                tempTrajStruct = trajStruct([trajStruct.task]==task);
                for i = 1:size(tempTrajStruct,2)
                    condObs = vertcat(tempTrajStruct(i).allSmoothFR.obs);
                    obs = vertcat(obs,condObs);
                    posture = tempTrajStruct(i).posture;
                    labels = vertcat(labels,ones(size(condObs,1),1)*posture);
                end
           end
           taskObsStruct(structInd).task = task;
           taskObsStruct(structInd).obs = obs;
           taskObsStruct(structInd).labels = labels;
           taskObsStruct(structInd).c = cvpartition(size(obs,1),'KFold',numFolds);
           structInd = structInd + 1;
        end

        %% Decode - train on one task, test on others
        result = struct('trainTask',[],'testTask',[],'pctCorrect',[],'avgPctCorrect',[]);
        %Set up result struct
        structInd = 1;
        numFolds = 5;
        for trainTask = 1:4
            for testTask = 1:4
                result(structInd).trainTask = trainTask;
                result(structInd).testTask = testTask;
                result(structInd).pctCorrect = zeros(1,numFolds);
                structInd = structInd + 1;
            end
        end

        for fold = 1:numFolds
            for trainTask = 1:4 
               trainC = taskObsStruct([taskObsStruct.task]==trainTask).c;
               trainIdx = training(trainC,fold);
               obs = taskObsStruct([taskObsStruct.task]==trainTask).obs;
               labels = taskObsStruct([taskObsStruct.task]==trainTask).labels;
               [LDAModel,numClasses] = LDATrainModel(obs(trainIdx,:),labels(trainIdx,1));
               %[NBModel,numClasses] = NBTrainModel(obs(trainIdx,:),labels(trainIdx,1));
               %Test model on testTask
               for testTask = 1:4
                    testC = taskObsStruct([taskObsStruct.task]==testTask).c;
                    testIdx = test(testC,fold);
                    obs = taskObsStruct([taskObsStruct.task]==testTask).obs;
                    labels = taskObsStruct([taskObsStruct.task]==testTask).labels;
                    [predictedLabels,foldPctCorrect] = LDAClassify(obs,labels,numClasses,LDAModel);
                    %[predictedLabels,foldPctCorrect] = NBClassify(obs,labels,numClasses,NBModel);
                    result([result.trainTask]==trainTask & [result.testTask]==testTask).pctCorrect(fold) = foldPctCorrect;
               end
           end
        end

        for i = 1:size(result,2)
           result(i).avgPctCorrect = mean(result(i).pctCorrect);     
        end

        resultStruct(resultStructInd).dataset = dataset;
        resultStruct(resultStructInd).result = result;
        resultStructInd = resultStructInd + 1;
    end
    
%% Plot performance
    figure; hold on
    numDatasets = size(resultStruct,2);
    plotInd = 1;
    for trainTask = 1:3
        for testTask = 1:3
            allPctCorrect = zeros(1,numDatasets);
            for i = 1:numDatasets
                result = resultStruct(i).result;
                allPctCorrect(i) = result([result.trainTask]==trainTask & [result.testTask]==testTask).avgPctCorrect;
            end
            plot(plotInd,allPctCorrect,'.','MarkerSize',10);
            plot(plotInd,mean(allPctCorrect),'.','MarkerSize',20);
            plotInd = plotInd + 1;
        end        
    end
    ax = gca;
    ax.YLim = [0,100];
    xLim = ax.XLim;
    plot(xLim,[100/3 100/3],'--','Color','r')
    
    ylabel('% Correct')
    fs = 14;
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)

    
    
%     xticklabelcell = cell(1,size(resultStruct,2));
%     for i = 1:size(resultStruct,2)
%        pctCorrect = resultStruct(i).pctCorrect; 
%        plot(i,pctCorrect,'.','MarkerSize',10,'Color','b')
%        plot(i,mean(pctCorrect),'.','MarkerSize',15,'Color','k')
%        trainTask = resultStruct(i).trainTask;
%        testTask = resultStruct(i).testTask;
%        xticklabelcell{1,i} = ['Train: ',taskIDs(trainTask).task,' Test: ',taskIDs(testTask).task];
%     end
%     
%     ax = gca;
%     xLim = ax.XLim;
%     plot(xLim,[100/3 100/3],'--','Color','r')
%     ax.YLim = [0,100];
%     ylabel('% Correct')
%     xticks([1:16])
%     xticklabels(xticklabelcell)
%     xtickangle(90)
%         fs = 14;
%     set(gca,'fontname','arial')
%     set(gca,'fontsize',fs)

