clear; clc; clf; close all;

%% Load Data 
    %3/17/2020
    date = '20200317';
    exclCh =  [3 12 16 20 31 73 85 92 94];
    binWidth = 1;
    [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
     N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200317(binWidth);
 
%% Collect States of Interest into Traj Field
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    numTrials = size(Data,2);
    for trial = 1:numTrials
        %Step 1
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
        [step1FR,step1FRTimestamps] = getStatesTraj(Data(trial),trialInclStates,'allChannelSpikeBins'); 
        Data(trial).Traj.step1FR = step1FR; Data(trial).Traj.step1FRTimestamps = step1FRTimestamps;
        [step1Marker,step1MarkerTimestamps] = getStatesTraj(Data(trial),trialInclStates,'marker'); 
        Data(trial).Traj.step1Marker = step1Marker; Data(trial).Traj.step1MarkerTimestamps = step1MarkerTimestamps;
    end
    clearvars trialInclStates
    
%% Keep only trials with length w/in 2 stdDevs of mean length
    [Data] = excludeLengths(Data,'step1FR');

%% Run GPFA On States of Interest; Save projection to Traj Field
    numDims = 12;
    binWidth = 25;
    D = struct('trialId',[],'spikes',[]);
    numTrials = size(Data,2);
    D = repmat(D,1,numTrials);
    for trial = 1:numTrials
       D(trial).trialId = trial;
       D(trial).spikes = Data(trial).Traj.step1FR';
    end
    % Run neuralTraj
    rmdir mat_results s
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
        GPFABinTimes = binWidth:binWidth:binWidth*GPFALength;
        Data(trial).Traj.step1GPFA = GPFAProj;
        Data(trial).Traj.step1GPFATimestamps = GPFABinTimes;
    end
    CSGPFAParams.Corth = Corth;
    clearvars result
    clearvars allSpikeBins GPFABinTimes GPFALength GPFAProj trialIdList TT Xorth Corth trialInd

%% Create trajStruct
    condFields = {{'posture','postureID'},{'target','targetData','target1ID'}};
    trajFields = {'step1FR','step1Marker','step1GPFA'};
    trajStruct = getTrajStruct(Data,condFields,trajFields);
 
%% Get "postural component" from trajectory IC's 
    postComp = struct('posture',[],'comp',[]);
    for posture = 1:5
        tempData = trajStruct([trajStruct.posture]==posture);
        numAvgTraj = size(tempData,2);
        ICs = zeros(numAvgTraj,numDims);
        for i = 1:numAvgTraj
            ICs(i,:) = tempData(i).avgStep1GPFA(1,:);
        end
        postComp(posture).posture = posture;
        postComp(posture).comp = mean(ICs);
    end
    
%% Get "target component" from target average trajectories
    targComp = struct('target',[],'comp',[]);
    for target = 1:8
        tempData = trajStruct([trajStruct.target]==target);
        numAvgTraj = size(tempData,2);
        trajLengths = zeros(1,numAvgTraj);
        for i = 1:numAvgTraj
           trajLengths(i) = size(tempData(i).avgStep1GPFA,1); 
        end
        allTargTraj = nan(max(trajLengths),numDims,numAvgTraj);
        for i = 1:numAvgTraj
            allTargTraj(1:trajLengths(i),:,i) = tempData(i).avgStep1GPFA;
        end
        avgTargTraj = nanmean(allTargTraj,3);
        avgTargTraj = avgTargTraj(1:round(mean(trajLengths)),:,:);
        targComp(target).target = target;
        targComp(target).comp = avgTargTraj;
    end
    
%% Create "modelled" trajectories
    modelStruct = struct('posture',[],'target',[],'traj',[]);
    structInd = 1;
    for posture = 1:5
        for target = 1:8
            condPostComp = postComp(posture).comp;
            condTargComp = targComp(target).comp;
            modelStruct(structInd).posture = posture;
            modelStruct(structInd).target = target;            
            modelStruct(structInd).traj = condTargComp + condPostComp;
            structInd = structInd + 1;
        end
    end
    

%% Use condition-averaged trajectories to perform posture, target, and PxT LDA
    avgTraj = vertcat(trajStruct.avgStep1GPFA);
    postureLabels = zeros(size(avgTraj,1),1);
    targetLabels = zeros(size(avgTraj,1),1);
    pxtLabels = zeros(size(avgTraj,1),1);
    labelInd = 1;
    for i = 1:size(trajStruct,2)
        postureLabel = trajStruct(i).posture;
        targetLabel = trajStruct(i).target;
        numPoints = size(trajStruct(i).avgStep1GPFA,1);
        postureLabels(labelInd:labelInd+numPoints-1,1) = postureLabel;
        targetLabels(labelInd:labelInd+numPoints-1,1) = targetLabel;
        pxtLabels(labelInd:labelInd+numPoints-1,1) = str2double([num2str(postureLabel),num2str(targetLabel)]);
        labelInd = labelInd + numPoints;
    end
    postureLDA = fisherLDA(avgTraj, postureLabels);
    [postureLDA,~] = qr(postureLDA);
    postureLDA = postureLDA(:,1:4);
    
    targetLDA = fisherLDA(avgTraj, targetLabels);
    [targetLDA,~] = qr(targetLDA);
    targetLDA = targetLDA(:,1:7);
    
    [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:3);
    
    pxtLDA = fisherLDA(avgTraj, pxtLabels);
    [pxtLDA,~] = qr(pxtLDA);
    pxtLDA = pxtLDA(:,1:numDims-1);

%% Add LDA Projections to Traj Struct and Model Struct
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgStep1PostureLDA = trajStruct(i).avgStep1GPFA*postureLDA;
       trajStruct(i).avgStep1TargetLDA = trajStruct(i).avgStep1GPFA*targetLDA;
       trajStruct(i).avgStep1PostTargOrth = trajStruct(i).avgStep1GPFA*postTargOrth;
       trajStruct(i).avgStep1PxTLDA = trajStruct(i).avgStep1GPFA*pxtLDA;
    end

    for i = 1:size(modelStruct,2)
       modelStruct(i).trajPostureLDA = modelStruct(i).traj*postureLDA;
       modelStruct(i).trajTargetLDA = modelStruct(i).traj*targetLDA;
       modelStruct(i).trajPostTargOrth = modelStruct(i).traj*postTargOrth;
       modelStruct(i).trajPxTLDA = modelStruct(i).traj*pxtLDA;
    end
    
%% Visualize Model and Actual Trajectories 
    load('C:\Users\pmari\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    dirStr = 'C:\Users\pmari\Documents\Posture\Figures\Decomposition Analysis\Earl20200317_Decomposition Analysis\';
    
    %Top 3 GPFA - Actual
    figure
    xDim = 1; yDim =2; zDim = 3;
    for posture = 1:5
        for target = 1:8
            traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgStep1GPFA;
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['GPFA ',num2str(xDim)]); ylabel(['GPFA ',num2str(yDim)]); zlabel(['GPFA ',num2str(zDim)])
    grid on
    titleStr = ['Actual - Top 3 GPFA']; title(titleStr);
    saveas(gcf,[dirStr,titleStr,'.fig'])
    
    %Top 3 GPFA - Model
    figure
    xDim = 1; yDim =2; zDim = 3;
    for posture = 1:5
        for target = 1:8
            traj = modelStruct(find([modelStruct.posture]==posture & [modelStruct.target]==target)).traj;
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['GPFA ',num2str(xDim)]); ylabel(['GPFA ',num2str(yDim)]); zlabel(['GPFA ',num2str(zDim)]);
    grid on
    titleStr = ['Model - Top 3 GPFA']; title(titleStr);
    saveas(gcf,[dirStr,titleStr,'.fig'])
%     %Posture/Target Projection - Actual
%     figure
%     xDim = 1; yDim = 1; zDim = 2;
%     for posture = 1:5
%         for target = 1:8
%             postureProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgStep1PostureLDA;
%             targetProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgStep1TargetLDA;
%             plot3(postureProj(:,xDim),targetProj(:,yDim),targetProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
%             hold on
%         end
%     end
%     xlabel(['Posture LDA ',num2str(xDim)]); ylabel(['Target LDA ',num2str(yDim)]); zlabel(['Target LDA ',num2str(zDim)])
%     grid on
%     titleStr = ['Actual']; title(titleStr);
%     
%     %Posture/Target Projection - Model
%     figure
%     xDim = 1; yDim = 1; zDim = 2;
%     for posture = 1:5
%         for target = 1:8
%             postureProj = modelStruct(find([modelStruct.posture]==posture & [modelStruct.target]==target)).trajPostureLDA;
%             targetProj = modelStruct(find([modelStruct.posture]==posture & [modelStruct.target]==target)).trajTargetLDA;
%             plot3(postureProj(:,xDim),targetProj(:,yDim),targetProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
%             hold on
%         end
%     end
%     xlabel(['Posture LDA ',num2str(xDim)]); ylabel(['Target LDA ',num2str(yDim)]); zlabel(['Target LDA ',num2str(zDim)])
%     grid on
%     titleStr = ['Model']; title(titleStr);
    
    %Posture/Target Orth Projection - Actual
    figure
    xDim = 1; yDim = 2; zDim = 3;
    for posture = 1:5
        for target = 1:8
            postTargOrthProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgStep1PostTargOrth;
            plot3(postTargOrthProj(:,xDim),postTargOrthProj(:,yDim),postTargOrthProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['Posture LDA 1']); ylabel(['Target LDA 1']); zlabel(['Target LDA 2'])
    grid on
    titleStr = ['Actual - Posture-arget']; title(titleStr);
    saveas(gcf,[dirStr,titleStr,'.fig'])
    
    %Posture/Target Orth Projection - Model
    figure
    xDim = 1; yDim = 2; zDim = 3;
    for posture = 1:5
        for target = 1:8
            postTargOrthProj = modelStruct(find([modelStruct.posture]==posture & [modelStruct.target]==target)).trajPostTargOrth;
            plot3(postTargOrthProj(:,xDim),postTargOrthProj(:,yDim),postTargOrthProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['Posture LDA 1']); ylabel(['Target LDA 1']); zlabel(['Target LDA 2'])
    grid on
    titleStr = ['Model - Posture-Target']; title(titleStr);
    saveas(gcf,[dirStr,titleStr,'.fig'])
    
    %PxT Projection - Actual
    figure
    xDim = 1; yDim = 2; zDim = 3;
    for posture = 1:5
        for target = 1:8
            pxtProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgStep1PxTLDA;
            plot3(pxtProj(:,xDim),pxtProj(:,yDim),pxtProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['PxT LDA ',num2str(xDim)]); ylabel(['PxT LDA ',num2str(yDim)]); zlabel(['PxT LDA ',num2str(zDim)])
    grid on
    titleStr = ['Actual - PxT']; title(titleStr);
    saveas(gcf,[dirStr,titleStr,'.fig'])
    
    %PxT Projection - Model
    figure
    xDim = 1; yDim = 2; zDim = 3;
    for posture = 1:5
        for target = 1:8
            pxtProj = modelStruct(find([modelStruct.posture]==posture & [modelStruct.target]==target)).trajPxTLDA;
            plot3(pxtProj(:,xDim),pxtProj(:,yDim),pxtProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['PxT LDA ',num2str(xDim)]); ylabel(['PxT LDA ',num2str(yDim)]); zlabel(['PxT LDA ',num2str(zDim)])
    grid on
    titleStr = ['Model - PxT']; title(titleStr);
    saveas(gcf,[dirStr,titleStr,'.fig'])
    
%% Assess Model Error 