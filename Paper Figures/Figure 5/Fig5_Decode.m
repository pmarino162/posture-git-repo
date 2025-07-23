clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 5';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Parameters
   numPCsToKeep = 10;     %Num PCs to project data into before analysis
   numFolds = 5;
   
%% Setup resultStruct
    resultStruct = struct('dataset',[],'result',[]);
    resultStructInd = 1;
    
%% Load and preprocess data
    %Load data
    dataset = 'E20200314';
    [Data,zScoreParams] = loadData(dataset);

    %Get trajStruct
    [Data] = removeShortBCIandIsoTrials(Data,dataset);
    [~,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);  
    [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 
    % Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(numPostures);          
    %For each timecourse, get one point
     for i = 1:size(trajStruct,2)
        for j = 1:size(trajStruct(i).allZSmoothFR,2)
           trajStruct(i).allZSmoothFR(j).obs = mean(trajStruct(i).allZSmoothFR(j).traj); 
        end
     end
    %Reduce dimensionality of observations
    allObs = [];
    for i = 1:size(trajStruct,2)
        for j = 1:size(trajStruct(i).allZSmoothFR,2)
            obs = trajStruct(i).allZSmoothFR(j).obs;
            allObs = vertcat(allObs,obs);
        end
    end
    [coeff,score,latent,tsquared,explained,mu] = pca(allObs); 
    for i = 1:size(trajStruct,2)
        for j = 1:size(trajStruct(i).allZSmoothFR,2)
            trajStruct(i).allZSmoothFR(j).obs = (trajStruct(i).allZSmoothFR(j).obs-mu)*coeff(:,1:numPCsToKeep);
        end
    end

%% Create taskObsStruct, balancing observations by target within each task x posture        
    %Get minimum number of trials to any target direction within each task x posture
    minNumObsStruct = struct('task',[],'posture',[],'minNumObs',[]);
    structInd = 1;
    for task = 1:3
        for posture = 1:3
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            tempMinNumCondObs = inf;
            for i = 1:size(tempTrajStruct,2)
                numCondObs = size(vertcat(tempTrajStruct(i).allZSmoothFR.obs),1);
                if numCondObs < tempMinNumCondObs
                    tempMinNumCondObs = numCondObs;
                end
            end
            minNumObsStruct(structInd).task = task;
            minNumObsStruct(structInd).posture = posture;
            minNumObsStruct(structInd).minNumObs = tempMinNumCondObs;
            structInd = structInd + 1;
        end
    end

    %Get max possible datsaset sizes for each task, choose smallest
    taskNumTargets = [8,3,2];
    maxPossibleSizes = zeros(1,3);
    for task = 1:3
       for posture = 1:3
          numObsToUsePerTarget = minNumObsStruct([minNumObsStruct.task]==task & [minNumObsStruct.posture]==posture).minNumObs;
          maxPossibleSizes(task) = maxPossibleSizes(task) + numObsToUsePerTarget*taskNumTargets(task);
       end
    end
    datasetSize = min(maxPossibleSizes);
        
    %Build taskObsStruct
    %For each task, subsample a dataset that is target-balanced for each posture. Then subsample randomly from this to datasetSize.
    %For 'all tasks' case, subsample training set for each fold equally from each task
    taskObsStruct = struct('task',[],'obs',[],'labels',[],'c',[]);
    structInd = 1;
    for task = 1:3
        obs = []; labels = [];
        for posture = 1:3
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            numObsToUsePerTarget = minNumObsStruct([minNumObsStruct.task]==task & [minNumObsStruct.posture]==posture).minNumObs;
            for i = 1:size(tempTrajStruct,2)
                condObs = vertcat(tempTrajStruct(i).allZSmoothFR.obs);
                obs = vertcat(obs,datasample(condObs,numObsToUsePerTarget,'Replace',false));
                posture = tempTrajStruct(i).posture;
                labels = vertcat(labels,ones(numObsToUsePerTarget,1)*posture);
            end  
        end
        [obs,idx] = datasample(obs,datasetSize,'Replace',false);
        labels = labels(idx,:);
        taskObsStruct(structInd).task = task;
        taskObsStruct(structInd).obs = obs;
        taskObsStruct(structInd).labels = labels;
        taskObsStruct(structInd).c = cvpartition(size(obs,1),'KFold',numFolds);
        structInd = structInd + 1;
    end
 
%% Decode
    %Set up result struct
    result = struct('trainTask',[],'testTask',[],'pctCorrect',[],'avgPctCorrect',[]);
    structInd = 1;
    for trainTask = 1:4
        for testTask = 1:3
            result(structInd).trainTask = trainTask;
            result(structInd).testTask = testTask;
            result(structInd).pctCorrect = zeros(1,numFolds);
            result(structInd).avgPctCorrect = 0;
            result(structInd).SEM = 0;
            structInd = structInd + 1;
        end
    end

    %Do within and across tasks cases
    for fold = 1:numFolds
        for trainTask = 1:3
           trainC = taskObsStruct([taskObsStruct.task]==trainTask).c;
           trainIdx = training(trainC,fold);
           obs = taskObsStruct([taskObsStruct.task]==trainTask).obs;
           labels = taskObsStruct([taskObsStruct.task]==trainTask).labels;
           [LDAModel,numClasses] = LDATrainModel(obs(trainIdx,:),labels(trainIdx,1));
           for testTask = 1:3
                testC = taskObsStruct([taskObsStruct.task]==testTask).c;
                testIdx = test(testC,fold);
                obs = taskObsStruct([taskObsStruct.task]==testTask).obs;
                labels = taskObsStruct([taskObsStruct.task]==testTask).labels;
                [predictedLabels,foldPctCorrect] = LDAClassify(obs(testIdx,:),labels(testIdx,1),numClasses,LDAModel);
                result([result.trainTask]==trainTask & [result.testTask]==testTask).pctCorrect(fold) = foldPctCorrect;
           end
       end
    end

    %Do 'train on all, test on all' case
    trainTask = 4;
    for fold = 1:numFolds
        %Create training dataset by subsampling from current training fold of each task
        numPerTask = floor(datasetSize/3);
        obs = NaN(numPerTask*3,numPCsToKeep);
        labels = NaN(numPerTask*3,1);
        for task = 1:3
           trainC = taskObsStruct([taskObsStruct.task]==task).c;
           trainIdx = training(trainC,fold);
           allObs = taskObsStruct([taskObsStruct.task]==task).obs(trainIdx,:);
           allLabels = taskObsStruct([taskObsStruct.task]==task).labels(trainIdx,1);
           [obs(numPerTask*(task-1)+1:numPerTask*task,:),idx] = datasample(allObs,numPerTask,'Replace',false);
           labels(numPerTask*(task-1)+1:numPerTask*task,:) = allLabels(idx,:);
        end
       [LDAModel,numClasses] = LDATrainModel(obs,labels);
       %Test model on testTask
       for testTask = 1:3
            testC = taskObsStruct([taskObsStruct.task]==testTask).c;
            testIdx = test(testC,fold);
            obs = taskObsStruct([taskObsStruct.task]==testTask).obs;
            labels = taskObsStruct([taskObsStruct.task]==testTask).labels;
            [predictedLabels,foldPctCorrect] = LDAClassify(obs(testIdx,:),labels(testIdx,1),numClasses,LDAModel);
            result([result.trainTask]==trainTask & [result.testTask]==testTask).pctCorrect(fold) = foldPctCorrect;
       end
    end
    
    %Get stats
    for i = 1:size(result,2)
       result(i).avgPctCorrect = mean(result(i).pctCorrect);     
    end

%% Plot performance - bar
withinTask = [];
acrossTask = [];
trainAllTestAll = [];

for i = 1:size(result,2)
    pctCorrect = result(i).pctCorrect;
    trainTask = result(i).trainTask;
    testTask = result(i).testTask;
        if trainTask < 4
           if testTask == trainTask
               withinTask = [withinTask,pctCorrect];
           else
               acrossTask = [acrossTask,pctCorrect];
           end
        else      
           trainAllTestAll = [trainAllTestAll,pctCorrect];
        end
end

figure; hold on
for barNum = 1:3
   if barNum == 1
       tempData = withinTask;
   elseif barNum == 2
       tempData = trainAllTestAll;
   elseif barNum == 3
       tempData = acrossTask;
   end
        tempMean = mean(tempData);
        SEM =  std(tempData)/sqrt(length(tempData));
        bar(barNum,mean(tempData),'w');
        errorbar(barNum,tempMean,SEM,'k','LineWidth',3);
end

ax = gca;
    yticks([0 50 100])
    ylabel('Posture Classification Accuracy (%)');
    fs = 14;
    set(gca,'fontname','arial'); set(gca,'fontsize',fs);
    ax.YLim = [0,100];
    xLim = ax.XLim;
    plot(xLim,[100/3 100/3],'--','Color','r')
    ax.TickDir='out';