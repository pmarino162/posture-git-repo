clear; clc; clf; close all;

%% Load, Preprocess, and Label Data    
    [Data] = loadEarlData20211004;
    Data = Data([Data.trialStatus]==1);
    
%% Exlude trials that are too long
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};

%% Run GPFA On States of Interest; Save projection to Spikes Field
    %Info for getStatesTraj    
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Center Target','Step 1'};
        trialInclStates(1).inclOccurrence = {'last','last'};
        trialInclStates(1).addTimeToBeginning = {-100,0};
        trialInclStates(1).addTimeToEnd = {0,100};
        

    %Info for GPFA
    numDims = 17;
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
    condFields = {{'task','conditionData','taskID'},{'posture','conditionData','postureID'},{'target','targetData','targetID'}};
    trajFields = {'GPFA','marker'};

    %Execution
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {0};
        trialInclStates(1).addTimeToEnd = {0};
        
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
            

%% BC Target & Posture LDA
    task = 1;
    tempTrajStruct = trajStruct([trajStruct.task]==task);
    postureList = unique([tempTrajStruct.posture]);
    targetList = unique([tempTrajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    allTraj = []; allPostureLabels = []; allTargetLabels = [];
    for i = 1:size(tempTrajStruct,2)
       target = tempTrajStruct(i).target;
       posture = tempTrajStruct(i).posture;
       numTraj = size(tempTrajStruct(i).allGPFA,2);
       for j = 1:numTraj
            %Trajectory
            traj = tempTrajStruct(i).allGPFA(j).traj;
            numSteps = size(traj,1);
            allTraj = vertcat(allTraj,traj);
            %Labels
            trajTargetLabel = ones(numSteps,1).*target;
            allTargetLabels = vertcat(allTargetLabels,trajTargetLabel);
            trajPostureLabel = ones(numSteps,1).*posture;
            allPostureLabels = vertcat(allPostureLabels,trajPostureLabel);
        end
    end
    
    
    %All Posture and Target
    postureLDA = fisherLDA(allTraj, allPostureLabels);
    [postureLDA,~] = qr(postureLDA);
    BCPostureLDA = postureLDA(:,1:numPostures-1);
    
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    BCTargetLDA = targetLDA(:,1:numTargets-1);
    
    [postTargOrth,~] = qr([BCPostureLDA(:,1:3),BCTargetLDA(:,1:2)]);
    BCPostTargOrth = postTargOrth(:,1:5);
    

%% Add LDA Projections to Traj Struct and Model Struct
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgBCPostureLDA = trajStruct(i).avgGPFA.traj*BCPostureLDA;
       trajStruct(i).avgBCTargetLDA = trajStruct(i).avgGPFA.traj*BCTargetLDA;
       trajStruct(i).avgBCPostTargOrth = trajStruct(i).avgGPFA.traj*BCPostTargOrth;      
    end

%% Load Colormaps
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customBlue.mat');
    ecmap = customBlue;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    
%% 3D Views of BCI Trajectories
	task  = 1;     
    tempTrajStruct = trajStruct([trajStruct.task]==task);
    postureList = unique([tempTrajStruct.posture]);
    targetList = unique([tempTrajStruct.target]);
    
    %Posture Space
    figure
    xDim = 1; yDim = 2; zDim = 3;
    for posture = postureList
        for target = targetList
            postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgBCPostTargOrth;
            plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
            plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
        end
    end
    xlabel(['Posture LDA 1']); ylabel(['Posture LDA 2']); zlabel(['Posture LDA 3'])
    grid on
    axis equal

%     posture = 2;
%     tempTrajStruct = delayTrajStruct;
%     for target = targetList
%         postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgBCPostTargOrth;
%         plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',scmap(3,:),'LineWidth',2); 
%         hold on
%         plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
%         plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
%         plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
%      end
    
    %% Inidividual Posture + Target Plane
    for postureDim = 1:3
        figure
        xDim = postureDim; yDim = 4; zDim = 5;
        for posture = postureList
            switch posture
                case 1
                    cmap = scmap(1,:);
                case 3
                    cmap = scmap(5,:);
                case 4
                    cmap = ecmap(1,:);
                case 5
                    cmap = ecmap(5,:);
            end
            for target = targetList
                postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgBCPostTargOrth;
                plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(1,:),'LineWidth',2); 
                hold on
                plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
            end
        end
        xlabel(['Posture LDA ',num2str(postureDim)]); ylabel(['Target LDA 1']); zlabel(['Target LDA 2'])
        grid on
        axis equal
    end