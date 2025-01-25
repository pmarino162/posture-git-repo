clear; clc; clf; close all;

%% Load Data 
    %Earl - 3/17/2020
    date = '20200317';
    exclCh =  [3 12 16 20 31 73 85 92 94];
    binWidth = 1;
    [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
     N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200317(binWidth);

%% Keep only trials with length w/in 2 stdDevs of mean length
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    trialInclStates(1).inclStates = {'Step 1'}; trialInclStates(1).inclOccurrence = {'last'};
    
    [Data] = excludeLengths(Data,trialInclStates);    
    
%% Run GPFA On States of Interest; Save projection to Traj Field
    numDims = 10;
    binWidth = 45;
    D = struct('trialId',[],'spikes',[]);
    numTrials = size(Data,2);
    D = repmat(D,1,numTrials);
    for trial = 1:numTrials
       D(trial).trialId = trial;
       D(trial).spikes = double(Data(trial).spikes.allChannelSpikeBins');
    end
    % Run neuralTraj
%     rmdir mat_results s
    result = neuralTraj(001,D,'binWidth',binWidth,'xDim',numDims); 
    CSGPFAParams = result;
    %[estParams, seqTrain, seqTest] = postprocess(result, varargin);
    clearvars D
    % Orthonormalize GPFA results and add to Data Struct
    trialIdList = [result.seqTrain.trialId];
    for trial = 1:numTrials
        trialInd = find(trialIdList == trial);
        [Xorth, Corth, TT]= orthogonalize(result.seqTrain(trialInd).xsm,result.estParams.C);
        GPFAProj = Xorth';
        GPFALength = result.seqTrain(trialInd).T;
        binStartTime = 1;
        binEndTime = binStartTime+binWidth*(GPFALength-1);
        GPFABinTimes = binStartTime:binWidth:binEndTime;
        Data(trial).spikes.GPFAProj = GPFAProj;
        Data(trial).spikes.GPFABinTimes = GPFABinTimes;
    end
    CSGPFAParams.Corth = Corth;
    clearvars result
    clearvars allSpikeBins GPFABinTimes GPFALength GPFAProj trialIdList TT Xorth Corth trialInd
    rmdir mat_results s
    
%% Collect States of Interest into Traj Field
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    trialInclStates(1).inclStates = {'Step 1'}; trialInclStates(1).inclOccurrence = {'last'};
    trialInclStates(1).addTimeToBeginning = {-500};
    
    numTrials = size(Data,2);
    for trial = 1:numTrials
        [GPFA,GPFATimestamps] = getStatesTraj(Data(trial),trialInclStates,'GPFAProj');
        Data(trial).Traj.GPFA = GPFA; Data(trial).Traj.GPFATimestamps = GPFATimestamps;
    end
    
    

%% Create trajStruct
    condFields = {{'posture','postureData','postureID'},{'target','targetData','target1ID'}};
    trajFields = {'GPFA'};
    trajStruct = getTrajStruct(Data,condFields,trajFields,'addTimeToBeginning',-100);
    trajStruct = getTrajStruct(Data,condFields,trajFields);
    
    

%% Use all data to perform posture, target, and PxT LDA
allTraj = []; allPostureLabels = []; allTargetLabels = []; allPxTLabels = [];

for i = 1:size(trajStruct,2)
   target = trajStruct(i).target;
   posture = trajStruct(i).posture;
   numTraj = size(trajStruct(i).allGPFA,2);
   for j = 1:numTraj
        %Trajectory
        traj = trajStruct(i).allGPFA(j).traj;
        numSteps = size(traj,1);
        allTraj = vertcat(allTraj,traj);
        %Labels
        trajTargetLabel = ones(numSteps,1).*target;
        allTargetLabels = vertcat(allTargetLabels,trajTargetLabel);
        trajPostureLabel = ones(numSteps,1).*posture;
        allPostureLabels = vertcat(allPostureLabels,trajPostureLabel);
        allPxTLabels = vertcat(allPxTLabels,ones(numSteps,1).*str2double([num2str(posture),num2str(target)]));
    end
end
    postureLDA = fisherLDA(allTraj, allPostureLabels);
    [postureLDA,~] = qr(postureLDA);
    postureLDA = postureLDA(:,1:4);
    
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    targetLDA = targetLDA(:,1:7);
    
    [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:3);
    
    pxtLDA = fisherLDA(allTraj, allPxTLabels);
    [pxtLDA,~] = qr(pxtLDA);
    pxtLDA = pxtLDA(:,1:numDims-1);
    
% %% Use condition-averaged trajectories to perform posture, target, and PxT LDA
%     avgTraj = vertcat(trajStruct.avgGPFA);
%     avgTraj = vertcat(avgTraj.traj);
%     postureLabels = zeros(size(avgTraj,1),1);
%     targetLabels = zeros(size(avgTraj,1),1);
%     pxtLabels = zeros(size(avgTraj,1),1);
%     labelInd = 1;
%     for i = 1:size(trajStruct,2)
%         postureLabel = trajStruct(i).posture;
%         targetLabel = trajStruct(i).target;
%         numPoints = size(trajStruct(i).avgGPFA.traj,1);
%         postureLabels(labelInd:labelInd+numPoints-1,1) = postureLabel;
%         targetLabels(labelInd:labelInd+numPoints-1,1) = targetLabel;
%         pxtLabels(labelInd:labelInd+numPoints-1,1) = str2double([num2str(postureLabel),num2str(targetLabel)]);
%         labelInd = labelInd + numPoints;
%     end
%     postureLDA = fisherLDA(avgTraj, postureLabels);
%     [postureLDA,~] = qr(postureLDA);
%     postureLDA = postureLDA(:,1:4);
%     
%     targetLDA = fisherLDA(avgTraj, targetLabels);
%     [targetLDA,~] = qr(targetLDA);
%     targetLDA = targetLDA(:,1:7);
%     
%     [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
%     postTargOrth = postTargOrth(:,1:3);
%     
%     pxtLDA = fisherLDA(avgTraj, pxtLabels);
%     [pxtLDA,~] = qr(pxtLDA);
%     pxtLDA = pxtLDA(:,1:numDims-1);

%% Add LDA Projections to Traj Struct and Model Struct
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgPostureLDA = trajStruct(i).avgGPFA.traj*postureLDA;
       trajStruct(i).avgTargetLDA = trajStruct(i).avgGPFA.traj*targetLDA;
       trajStruct(i).avgPostTargOrth = trajStruct(i).avgGPFA.traj*postTargOrth;
       trajStruct(i).avgPxTLDA = trajStruct(i).avgGPFA.traj*pxtLDA;
    end

    
%% Visualize Model and Actual Trajectories 
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Decomposition Analysis\Earl20200317_Decomposition Analysis\';
    
    %Top 3 GPFA - Actual
    figure
    xDim = 1; yDim =2; zDim = 3;
    for posture = 1:5
        for target = 1:8
            traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgGPFA.traj;
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['GPFA ',num2str(xDim)]); ylabel(['GPFA ',num2str(yDim)]); zlabel(['GPFA ',num2str(zDim)])
    grid on
    titleStr = ['Actual - Top 3 GPFA']; title(titleStr);
%     saveas(gcf,[dirStr,titleStr,'.fig'])

    %Posture/Target Projection - Actual
    figure
    xDim = 1; yDim = 1; zDim = 2;
    for posture = 1:5
        for target = 1:8
            postureProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPostureLDA;
            targetProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgTargetLDA;
            plot3(postureProj(:,xDim),targetProj(:,yDim),targetProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['Posture LDA ',num2str(xDim)]); ylabel(['Target LDA ',num2str(yDim)]); zlabel(['Target LDA ',num2str(zDim)])
    grid on
    titleStr = ['Actual Posture-Target']; title(titleStr);
   
    %Posture/Target Orth Projection - Actual
    figure
    xDim = 1; yDim = 2; zDim = 3;
    for posture = 1:5
        for target = 1:8
            postTargOrthProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPostTargOrth;
            plot3(postTargOrthProj(:,xDim),postTargOrthProj(:,yDim),postTargOrthProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['Posture LDA 1']); ylabel(['Target LDA 1']); zlabel(['Target LDA 2'])
    grid on
    titleStr = ['Actual - Posture-Target Orth']; title(titleStr);
% %     saveas(gcf,[dirStr,titleStr,'.fig'])
%     

    %PxT Projection - Actual
    figure
    xDim = 1; yDim = 2; zDim = 3;
    for posture = 1:5
        for target = 1:8
            pxtProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPxTLDA;
            plot3(pxtProj(:,xDim),pxtProj(:,yDim),pxtProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['PxT LDA ',num2str(xDim)]); ylabel(['PxT LDA ',num2str(yDim)]); zlabel(['PxT LDA ',num2str(zDim)])
    grid on
    titleStr = ['Actual - PxT']; title(titleStr);
