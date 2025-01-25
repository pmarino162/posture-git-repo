clear; clc; clf; close all;

%% Load, Preprocess, and Label Data    
    [Data] = loadEarlData20210825;
    Data = Data([Data.trialStatus]==1);
    
%% Exlude trials that are too long
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
    trialInclStates(2).trialName = {'GridReaching'};
        trialInclStates(2).inclStates = {'No Cheating','Target Acquire'};
        trialInclStates(2).inclOccurrence = {'last','last'};
    [Data] = excludeLengths(Data,trialInclStates);

%% Run GPFA On States of Interest; Save projection to Spikes Field
    %Info for getStatesTraj    
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {-100};
        trialInclStates(1).addTimeToEnd = {100};
    trialInclStates(2).trialName = {'GridReaching'};
        trialInclStates(2).inclStates = {'No Cheating','Target Acquire'};
        trialInclStates(2).inclOccurrence = {'last','last'};
        trialInclStates(2).addTimeToBeginning = {-100,0,0};
        trialInclStates(2).addTimeToEnd = {0,100};

    %Info for GPFA
    numDims = 15;
    binWidth = 25;
    D = struct('trialId',[],'spikes',[]);
    numTrials = size(Data,2);
    D = repmat(D,1,numTrials);
    startTimes = zeros(1,numTrials);
    for trial = 1:numTrials
       D(trial).trialId = trial;
       [FR,FRTimestamps] = getStatesTraj(Data(trial),trialInclStates,'allChannelSpikeBins','timeRelToTrialStart',true);
       D(trial).spikes = FR';
       startTimes(trial) = FRTimestamps(1);
    end
    % Run neuralTraj
    result = neuralTraj(001,D,'binWidth',binWidth,'xDim',numDims); 
    CSGPFAParams = result;
    clearvars D
    % Orthonormalize GPFA results and add to Data Struct
    trialIdList = [result.seqTrain.trialId];
    for trial = 1:numTrials
        trialInd = find(trialIdList == trial);
        [Xorth, Corth, TT]= orthogonalize(result.seqTrain(trialInd).xsm,result.estParams.C);
        GPFAProj = Xorth';
        GPFALength = result.seqTrain(trialInd).T;
        binStartTime = startTimes(trial) + binWidth/2;
        binEndTime = binStartTime+binWidth*(GPFALength-1);
        GPFABinTimes = binStartTime:binWidth:binEndTime;
        Data(trial).spikes.GPFA = GPFAProj;
        Data(trial).spikes.GPFABinTimes = GPFABinTimes;
    end
    CSGPFAParams.Corth = Corth;
    clearvars result
    clearvars allSpikeBins GPFABinTimes GPFALength GPFAProj trialIdList TT Xorth Corth trialInd
    rmdir mat_results s
    
%% Create trajStruct
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'}};
    trajFields = {'GPFA','marker','absoluteForce'};
    
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {0};
        trialInclStates(1).addTimeToEnd = {0};
    trialInclStates(2).trialName = {'GridReaching'};
        trialInclStates(2).inclStates = {'No Cheating','Target Acquire'};
        trialInclStates(2).inclOccurrence = {'last','last','last'};
        trialInclStates(2).addTimeToBeginning = {0,0};
        trialInclStates(2).addTimeToEnd = {0,0};    
        
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
            
    
%% Use BC data to perform targetLDA
    BCTrajStruct = trajStruct([trajStruct.task]==1);
    allTraj = []; allTargetLabels = []; 
    for i = 1:size(BCTrajStruct,2)
       target = BCTrajStruct(i).target;
       numTraj = size(BCTrajStruct(i).allGPFA,2);
       for j = 1:numTraj
            %Trajectory
            traj = BCTrajStruct(i).allGPFA(j).traj;
            numSteps = size(traj,1);
            allTraj = vertcat(allTraj,traj);
            %Labels
            trajTargetLabel = ones(numSteps,1).*target;
            allTargetLabels = vertcat(allTargetLabels,trajTargetLabel);
        end
    end
    
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    BCTargetLDA = targetLDA(:,1:2);

    nullSpace = null(BCTargetLDA');
    nullProj = allTraj*nullSpace;
    nullSpacePC = pca(nullProj);
    BCTargetNullSpace = nullSpace*nullSpacePC;

    
%% Use Arm data to perform targetLDA
    HCTrajStruct = trajStruct([trajStruct.task]==2);
    allTraj = []; allTargetLabels = []; 
    for i = 1:size(HCTrajStruct,2)
       target = HCTrajStruct(i).target;
       numTraj = size(HCTrajStruct(i).allGPFA,2);
       for j = 1:numTraj
            %Trajectory
            traj = HCTrajStruct(i).allGPFA(j).traj;
            numSteps = size(traj,1);
            allTraj = vertcat(allTraj,traj);
            %Labels
            trajTargetLabel = ones(numSteps,1).*target;
            allTargetLabels = vertcat(allTargetLabels,trajTargetLabel);
        end
    end
    
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    HCTargetLDA = targetLDA(:,1:2);

    nullSpace = null(HCTargetLDA');
    nullProj = allTraj*nullSpace;
    nullSpacePC = pca(nullProj);
    HCTargetNullSpace = nullSpace*nullSpacePC;
    
%% Add LDA Projections to Traj Struct and Model Struct
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgBCTargetLDA = trajStruct(i).avgGPFA.traj*BCTargetLDA;
       trajStruct(i).avgHCTargetLDA = trajStruct(i).avgGPFA.traj*HCTargetLDA;
%        trajStruct(i).avgNull = trajStruct(i).avgGPFA.traj*nullSpace;
    end

%% Visualize Model and Actual Trajectories 
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    cmap = customRainbow;
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Decomposition Analysis\Earl20200317_Decomposition Analysis\';
    
%% Project BCI and Arm Data into BCI Target Plane
    for task = 1:2
        tempTrajStruct = trajStruct([trajStruct.task]==task);
        targetList = unique([tempTrajStruct.target]);
        
        figure 
        xDim = 1; yDim = 2;
        for target = targetList
            traj = tempTrajStruct([tempTrajStruct.target]==target).avgBCTargetLDA;
            plot(traj(:,xDim),traj(:,yDim),'Color',cmap(target,:),'LineWidth',2); 
            hold on
            plot(traj(1,xDim),traj(1,yDim),'.','Color',cmap(target,:),'MarkerSize',20);
        end
        xlabel('BC Target LDA 1')
        ylabel('BC Target LDA 1')
        
        switch task
            case 1
                title('BC Data')
            case 2
                title('HC Data')
        end
    end
    
%% Project BCI and Arm Data into Arm Target Plane 
    for task = 1:2
        tempTrajStruct = trajStruct([trajStruct.task]==task);
        targetList = unique([tempTrajStruct.target]);
        
        figure 
        xDim = 1; yDim = 2;
        for target = targetList
            traj = tempTrajStruct([tempTrajStruct.target]==target).avgHCTargetLDA;
            plot(traj(:,xDim),traj(:,yDim),'Color',cmap(target,:),'LineWidth',2); 
            hold on
            plot(traj(1,xDim),traj(1,yDim),'.','Color',cmap(target,:),'MarkerSize',20);
        end
        
        xlabel('HC Target LDA 1')
        ylabel('HC Target LDA 1')
        
        switch task
            case 1
                title('BC Data')
            case 2
                title('HC Data')
        end
        
    end
    

%% Get Forces During BCI 
    tempTrajStruct = trajStruct([trajStruct.task]==1);
    targetList = unique([tempTrajStruct.target]);
    figure
    
    for target = targetList
            traj = tempTrajStruct([tempTrajStruct.target]==target).avgAbsoluteForce.traj;
            time = tempTrajStruct([tempTrajStruct.target]==target).avgAbsoluteForce.timestamps;
            for dim = 1:3
                subplot(3,1,dim)
                plot(time,traj(:,dim),'Color',cmap(target,:),'LineWidth',2); 
                hold on
            end
            
    end
    
    for dim = 1:3
        subplot(3,1,dim)
        xlabel('time (ms)')
        switch dim
            case 1
                ylabel('F_X (N)')
            case 2
                ylabel('F_Y (N)')
            case 3
                ylabel('F_Z (N)')
        end
    end
    
    
    figure
    for target = targetList
        traj = tempTrajStruct([tempTrajStruct.target]==target).avgAbsoluteForce.traj;
        plot(traj(:,2),traj(:,3),'Color',cmap(target,:),'LineWidth',2); 
        hold on
    end
    axis equal
    xlabel('F_Y (N)')
    ylabel('F_Z (N)')
    
    
    
    figure
    for target = targetList
        traj = tempTrajStruct([tempTrajStruct.target]==target).avgAbsoluteForce.traj;
        plot(traj(:,2),traj(:,1),'Color',cmap(target,:),'LineWidth',2); 
        hold on
    end
    axis equal
    xlabel('F_Y (N)')
    ylabel('F_X (N)')
    
    figure
    for target = targetList
        traj = tempTrajStruct([tempTrajStruct.target]==target).avgAbsoluteForce.traj;
        plot3(traj(:,1),traj(:,2),traj(:,3),'Color',cmap(target,:),'LineWidth',2); 
        hold on
        plot3(traj(1,1),traj(1,2),traj(1,3),'.','Color',cmap(target,:),'MarkerSize',20); 
    end
    axis equal
    xlabel('F_X (N)')
    ylabel('F_Y (N)')
    zlabel('F_Z (N)')