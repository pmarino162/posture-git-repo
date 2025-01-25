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
    [Data] = loadEarlData20210901_20211210;
    
%% Get trajStruct
    binWidth = 25;
    trajFields = {'smoothFR','marker'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);

    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
    trialInclStates(2).trialName = {'DelayedCenterOut20210828'};
        trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
    
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
    minNumTimestamps = zeros(1,2);
    for task = 1:2
       tempTrajStruct = trajStruct([trajStruct.task]==task);
       numTimestamps = [];
       for i = 1:size(tempTrajStruct,2)
            numTimestamps(i) = length(tempTrajStruct(i).avgSmoothFR.timestamps);
       end
       [minNumTimestamps(task),~] = min(numTimestamps);
    end


%% Look at reaching position trajectories to understand time
    figure; hold on
    task = 2;
    for posture = 2
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj; 
            time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps;
            subplot(2,1,1)
                plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            subplot(2,1,2)
                plot(time,traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
        end
    end
    subplot(2,1,1)
        ylabel('x (mm)')
    subplot(2,1,2)
        xlabel('time (ms)')
        ylabel('y (mm)')
    
        
%% Form allTraj containing trial-averaged data for each condition
    allTraj = []; allPos = [];
    for i = 1:size(trajStruct,2)
        task = trajStruct(i).task;
        traj = trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps(task),:);
        allTraj = vertcat(allTraj,traj);
        pos = trajStruct(i).avgMarker.traj(1:minNumTimestamps(task),1);
        allPos = vertcat(allPos,pos);
    end

%% Do PCA on all data
    [PCA,score,latent,tsquared,explained,mu] = pca(allTraj); 
    
%% Do posture LDA on BCI
    taskTrajStruct = trajStruct([trajStruct.task]==1);
%     minNumTimestamps = minNumTimestamps(1);
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = [1 3 4 5]
       obsStruct(structInd).label = posture;
       allObs = [];
       tempTrajStruct = taskTrajStruct([taskTrajStruct.posture]==posture);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               traj = tempTrajStruct(i).allSmoothFR(j).traj;
               if size(traj,1) > minNumTimestamps(1)
                    traj = traj(1:minNumTimestamps(1),:);
               end
               allObs = vertcat(allObs,traj);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    [postureLDA] = doLDA(obsStruct);
    
%% Do linear regression to predict hand position from neural PCs
    allTraj = []; allPos = [];
    for i = 1:size(trajStruct,2)
        traj = trajStruct(i).avgSmoothFR(j).traj;
        pos = trajStruct(i).avgMarker(j).traj(:,1); 
        numTimestamps = size(pos,1);
        if numTimestamps > minNumTimestamps(task)
            traj = traj(1:minNumTimestamps(task),:);
            pos = pos(1:minNumTimestamps(task),:);
        end
            allTraj = vertcat(allTraj,traj);
            allPos = vertcat(allPos,pos);
    end

    X = allTraj; y = allPos;
    n = length(y);
    c = cvpartition(n,'HoldOut',0.3);
    idxTrain = training(c,1);
    idxTest = ~idxTrain;
    
    kList = [0.01,0.1,1,10,25,50,100];
    result = zeros(length(kList),2);
    kInd = 1;
    for k = kList
        b = ridge(y(idxTrain),X(idxTrain,:),k,0);
        yhat = b(1) + X(idxTest,:)*b(2:end); 
        error = sum((y(idxTest)-yhat).^2);
        result(kInd,1) = k;
        result(kInd,2) = error;
        kInd = kInd + 1;
    end
    
    figure
    plot(result(:,1),result(:,2))
    
    figure
    k = 50;
    b = ridge(y(idxTrain),X(idxTrain,:),k,0);
    yhat = b(1) + X(idxTest,:)*b(2:end);  
    scatter(y(idxTest),yhat)
    hold on
    plot(y(idxTest),y(idxTest))
    xlabel('Actual x Pos')
    ylabel('Predicted x Pos')
    hold off
    postureReg = b(2:end);
    
%% Orthonormalize combinations of axes
    [PPCOrth,~] = qr([postureLDA(:,1),PCA(:,1:2)]); PPCOrth = PPCOrth(:,1:3);
    [PRPCOrth,~] = qr([postureReg(:,1),PCA(:,1:2)]); PRPCOrth = PRPCOrth(:,1:3);
%% Add projections to trajStruct
    for i = 1:size(trajStruct,2)
        %Add average traces to trajStruct
        trajStruct(i).PCA.traj = trajStruct(i).avgSmoothFR.traj*PCA;
        trajStruct(i).PPCOrth.traj = trajStruct(i).avgSmoothFR.traj*PPCOrth;
        trajStruct(i).PRPCOrth.traj = trajStruct(i).avgSmoothFR.traj*PRPCOrth;
        trajStruct(i).postureLDA.traj = trajStruct(i).avgSmoothFR.traj*postureLDA;
        trajStruct(i).postureReg.traj = trajStruct(i).avgSmoothFR.traj*postureReg;
    end
    
%% Look at some trajectories
    figure
    hold on
    timePts = 1:10;
    for task = 2
       for posture = [2]
           for target = [1,5]
               if any([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target) 
                    traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).PRPCOrth.traj(timePts,:); 
%                     traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(timePts,:);
%                     traj = traj-CIMargTraj;
%                     traj = traj*PPCOrth;
                    
                    plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               end
           end
       end
    end
    grid on
    

%% Look at reaching position trajectories to understand time
    figure; hold on
    task = 2;
    for posture = [2]
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).postureReg.traj; 
            time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
            subplot(3,1,1)
                plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj; 
            time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps;
            subplot(3,1,2)
                plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            subplot(3,1,3)
                plot(time,traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
        end
    end
    subplot(3,1,2)
        ylabel('x (mm)')
    subplot(3,1,3)
        xlabel('time (ms)')
        ylabel('y (mm)')    

%% Plot marker pose vs time
    figure; hold on
    timeInd = 1;
    task = 2;
    timePts = 1:18;
    plotTimePts = timePts;
    numTimestamps = 18;
    for posture = postureList
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).PPCOrth.traj(timePts,:); 
%             time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(timePts,:);
            plot(plotTimePts,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
            plotTimePts = plotTimePts + numTimestamps;
        end
    end
    
%     figure; hold on
    timeInd = 1;
    task = 2;
    timePts = 1:18;
    plotTimePts = timePts;
    numTimestamps = 18;
    for posture = postureList
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj(timePts,:); 
%             time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(timePts,:);
            plot(plotTimePts,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
            plotTimePts = plotTimePts + numTimestamps;
        end
    end
    
    %% Plot marker pose vs time
    figure; hold on
    timeInd = 1;
    task = 2;
    timePts = 1:18;
    plotTimePts = timePts;
    numTimestamps = 18;
    for posture = postureList
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(timePts,:); 
            traj = b(1) + traj*b(2:end);  
            plot(plotTimePts,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
            plotTimePts = plotTimePts + numTimestamps;
        end
    end
    
%     figure; hold on
    timeInd = 1;
    task = 2;
    timePts = 1:18;
    plotTimePts = timePts;
    numTimestamps = 18;
    for posture = postureList
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj(timePts,:); 
%             time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(timePts,:);
            plot(plotTimePts,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
            plotTimePts = plotTimePts + numTimestamps;
        end
    end
    
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