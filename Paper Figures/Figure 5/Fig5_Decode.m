clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Error Reduction Figures';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat');
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\royg.mat');
    roygCmap = royg;
    
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
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        [~,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);  
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 

        taskIDs = struct('ID',[],'task','');
        taskIDs(1).ID = 1; taskIDs(1).task = 'BC';
        taskIDs(2).ID = 2; taskIDs(2).task = 'HC';
        taskIDs(3).ID = 3; taskIDs(3).task = 'Iso';
        taskIDs(4).ID = 4; taskIDs(4).task = 'All';
        
        %For each timecourse, get one point
         for i = 1:size(trajStruct,2)
            %Individual trials
            for j = 1:size(trajStruct(i).allZSmoothFR,2)
               trajStruct(i).allZSmoothFR(j).obs = mean(trajStruct(i).allZSmoothFR(j).traj); 
            end
            %Condition average
            trajStruct(i).avgZSmoothFR.condAvg = mean(trajStruct(i).avgZSmoothFR.traj);
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
                trajStruct(i).allZSmoothFR(j).obs = trajStruct(i).allZSmoothFR(j).obs*coeff(:,1:numManPCs);
            end
        end

        %Get posture list
        postureList = unique([trajStruct.posture]);
    
        %% Create taskObsStruct, balancing observations by target within each task x posture
       
        %Get minimum number of trials to any target direction within each
        %task x posture
        minNumObsStruct = struct('task',[],'posture',[],'minNumObs',[]);
        structInd = 1;
        for task = 1:3
            for posture = 1:3
                tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
                tempMinNumCondObs = 999999999;
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

        %Build taskObsStruct
        taskObsStruct = struct('task',[],'obs',[],'labels',[],'c',[]);
        structInd = 1;
        for task = 1:4
            obs = []; labels = [];
            if task < 4
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

            elseif task == 4
               minNumTaskObs = 999999999999;
               for tempTask = 1:3
                   numTaskObs = size(taskObsStruct(tempTask).obs,1);
                   if numTaskObs < minNumTaskObs
                       minNumTaskObs = numTaskObs;
                   end
               end
               for tempTask = 1:3
                    taskObs = taskObsStruct(tempTask).obs;
                    taskLabels = taskObsStruct(tempTask).labels;
                    [tempObs,idx] = datasample(taskObs,minNumTaskObs,'Replace',false);
                    tempLabels = taskLabels(idx);
                    obs = vertcat(obs,tempObs);
                    labels = vertcat(labels,tempLabels);
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

        %Visualize
        for fold = 1:numFolds
            f = figure;
            f.Position = [10 10 1000 500];
            for trainTask = 1:3
               %Fit LDA space
               trainC = taskObsStruct([taskObsStruct.task]==trainTask).c;
               trainIdx = training(trainC,fold);
               obs = taskObsStruct([taskObsStruct.task]==trainTask).obs;
               labels = taskObsStruct([taskObsStruct.task]==trainTask).labels;
               LDAproj = fisherLDA(obs(trainIdx,:), labels(trainIdx,1));
               [LDAproj,~] = qr(LDAproj);
               LDAproj = LDAproj(:,1:numClasses-1);
               for testTask = 1:3
                   testC = taskObsStruct([taskObsStruct.task]==testTask).c;
                   testIdx = test(testC,fold);
                   obs = taskObsStruct([taskObsStruct.task]==testTask).obs(testIdx,:);
                   labels = taskObsStruct([taskObsStruct.task]==testTask).labels(testIdx,1);
                   
                   %Plot in LDA space
                   subplot(3,3,(trainTask-1)*3+testTask)
                   hold on;
                        for posture = 1:3
                            postureObs = obs(labels==posture,:);
                            ldaProj = postureObs*LDAproj;
                            plot(ldaProj(:,1),ldaProj(:,2),'.','MarkerSize',10,'Color',pcmap(posture,:));
                        end
                   
                    
               end
           end
           sgtitle(['Fold ',num2str(fold)])
        end
        
        resultStruct(resultStructInd).dataset = dataset;
        resultStruct(resultStructInd).result = result;
        resultStructInd = resultStructInd + 1;
    end
    
%% Plot performance
    figure; hold on
    numDatasets = size(resultStruct,2);
    plotInd = 1;
    for trainTask = 1:4
        for testTask = 1:4
            allPctCorrect = zeros(1,numDatasets);
            for i = 1:numDatasets
                result = resultStruct(i).result;
                pctCorrect = result([result.trainTask]==trainTask & [result.testTask]==testTask).pctCorrect;
            end
            SEM = std(pctCorrect)/sqrt(length(pctCorrect));
            errorbar(plotInd,mean(pctCorrect),SEM,'k','LineWidth',3);
            plot(plotInd,mean(pctCorrect),'.','MarkerSize',20,'Color','k');
%             if testTask == 1
%                 
%             elseif task == 2
%                 
%             elseif task == 3
%                 
%             elseif task == 4
%                 
%             end
%             plot(plotInd,allPctCorrect,'.','MarkerSize',10,'Color','b');
%             plot(plotInd,mean(allPctCorrect),'.','MarkerSize',15,'Color','k');
       
            plotInd = plotInd + 1;
        end        
    end
    ax = gca;
    ax.XTick = 1:16;
    ax.XTickLabel = {};
    yticks([0 50 100])
    ylabel('Posture Classification Accuracy (%)');
    fs = 14;
    set(gca,'fontname','arial'); set(gca,'fontsize',fs);
    ax.YLim = [0,100];
    xLim = ax.XLim;
    plot(xLim,[100/3 100/3],'--','Color','r')
    ax.TickDir='out';


    

%% Plot classification performance - heatmap
resultMat = zeros(4,4);
for trainTask = 1:4
    for testTask = 1:4
        allPctCorrect = zeros(1,numDatasets);
        for i = 1:numDatasets
            result = resultStruct(i).result;
            allPctCorrect(i) = result([result.trainTask]==trainTask & [result.testTask]==testTask).avgPctCorrect;
        end
        resultMat(trainTask,testTask) = mean(allPctCorrect);
    end
end

    
    figure
    clims = [0,100];
    set(gca,'fontname','arial'); set(gca,'fontsize',fs);


%     tempJet = jet;
%     tempJet = tempJet(33:64,:);
%     colormap(flip(tempJet))
    colormap(royg)
    imagesc(resultMat,clims)
    xticks([1:4])
    xticklabels({'BCI','Reaching','Isometric Force','All'})
    xlabel('Testing Set')
    yticks([1:4])
    yticklabels({'BCI','Reaching','Isometric Force','All'})
    ylabel('Training Set')
    c = colorbar;
    
    c.Ticks = [0 25 33.33 50 75 100];
%     c.TickLabels = {'Chance (33.33)','50','75','100'};
    c.Label.String = 'Posture Classification Accuracy (%)';
    ax = gca;
    axis square
    ax.XTick = []; ax.YTick = [];
        ax.FontName = 'Arial';
    ax.FontSize = 14;
    
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

%% Plot performance - bar
withinTask = [];
acrossTask = [];
trainAllTestAll = [];

for i = 1:size(result,2)
    avgPctCorrect = result(i).avgPctCorrect;
    trainTask = result(i).trainTask;
    testTask = result(i).testTask;
    if testTask ~= 4
        if trainTask < 4
           if testTask == trainTask
               withinTask = [withinTask,avgPctCorrect];
           else
               acrossTask = [acrossTask,avgPctCorrect];
           end
        else      
           trainAllTestAll = [trainAllTestAll,avgPctCorrect];
        end
    end
end

figure; hold on
for barNum = 1:3
   if barNum == 1
       tempData = withinTask;
   elseif barNum == 2
       tempData = acrossTask;
   elseif barNum == 3
       tempData = trainAllTestAll;
   end
        SEM =  std(tempData)/sqrt(length(tempData));
        bar(barNum,mean(tempData),'w');
        errorbar(barNum,mean(tempData),SEM,'k','LineWidth',3);
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