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

%% Load Data
    [Data,~,~,~,~,~,~] = loadEarlData20200314_20211210;

%% Get trajStruct
    binWidth = 25;
    trajFields = {'smoothFR','marker'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);

    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    
    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
        trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',0},{'state','Step 2','first',0}};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(2).inclStates = {{'state','Reach','first',0},{'state','Hold','first',0}};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
    
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%% Get posture and target lists
    taskList = unique([trajStruct.task]);
    numTasks = size(taskList,2);
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);  
    
%% Get minimum number of timestamps in condition averages for each task
    minNumTimestamps = zeros(1,3);
    for task = 1:3
       tempTrajStruct = trajStruct([trajStruct.task]==task);
       numTimestamps = [];
       for i = 1:size(tempTrajStruct,2)
            numTimestamps(i) = length(tempTrajStruct(i).avgSmoothFR.timestamps);
       end
       [minNumTimestamps(task),~] = min(numTimestamps);
    end

%% Form allTraj containing trial-averaged data for each condition
    allTraj = []; 
    for i = 1:size(trajStruct,2)
        task = trajStruct(i).task;
        traj = trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps(task),:);
        allTraj = vertcat(allTraj,traj);
    end

%% Do PCA on all data
    [PCA,score,latent,tsquared,explained,mu] = pca(allTraj); 

%% Do task LDA
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for task = taskList
       obsStruct(structInd).label = task;
       allObs = [];
       tempTrajStruct = trajStruct([trajStruct.task]==task);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               traj = tempTrajStruct(i).allSmoothFR(j).traj;
               if size(traj,1) > minNumTimestamps
                    traj = traj(1:minNumTimestamps,:);
               end
               allObs = vertcat(allObs,traj);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    [taskLDA] = doLDA(obsStruct);
    
%% Add projections to trajStruct
    for i = 1:size(trajStruct,2)
        %Add average traces to trajStruct
        trajStruct(i).PCA.traj = trajStruct(i).avgSmoothFR.traj*PCA;
        trajStruct(i).taskLDA.traj = trajStruct(i).avgSmoothFR.traj*taskLDA;
    end
    
    
%% Visualize data in PC space
   figure
    hold on
    for task = taskList
       taskTrajStruct = trajStruct([trajStruct.task]==task);
       taskTargetList = [taskTrajStruct.target];
       for posture = [5]
           
           for target = taskTargetList
               if any([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target) 
                    timePts = 1:10;%minNumTimestamps(task);
                    traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).PCA.traj(timePts,:); 
                    plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               end
           end
       end
    end
    grid on
    axis equal

%% Visualize data in task LDA space
   figure
    hold on
    for task = taskList
       taskTrajStruct = trajStruct([trajStruct.task]==task);
       taskTargetList = [taskTrajStruct.target];
       for posture = 3
           for target = taskTargetList
               if any([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target) 
                    timePts = 1:10;%minNumTimestamps(task);
                    traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).taskLDA.traj(timePts,:); 
                    plot(traj(:,1),traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                    plot(traj(1,1),traj(1,2),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot(traj(end,1),traj(end,2),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               end
           end
       end
    end
    grid on
    axis equal

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
 