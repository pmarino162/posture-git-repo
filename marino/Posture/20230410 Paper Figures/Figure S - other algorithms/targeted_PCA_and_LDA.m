clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20230203\orthographic dpca pop traj';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Load dataset; get trajStruct
    dataset = 'E20200316';
    [Data,zScoreParams] = loadData(dataset);

    %Get trajStruct
    binWidth = 25;
    kernelStdDev = 25;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    switch dataset
        %BCI
        case {'E20200316','E20200317','E20200318'}
            trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
            condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
        case {'N20171215','N20180221'}
            trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
        case {'R20201020','R20201021'}
            trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
        case 'E20210901'
            taskID = [Data.conditionData]; taskID = [taskID.taskID];
            Data = Data(taskID==1);
            trialInclStates(1).trialName = {'BCI Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
        %Iso
        case 'E20200116'
            trialInclStates(1).trialName = {'IsometricForce_1D'};   
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',500}};
        %Reaching
        case 'E20210706'
            trialInclStates(1).trialName = {'GridReaching'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
        case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
            trialInclStates(1).trialName = {'Nigel Dissociation'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
        case{'R20200221','R20200222'}
            trialInclStates(1).trialName = {'Rocky Dissociation'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

    %Keep postures 1-4 only for earl
    switch dataset
        case {'E20210706'}
                trajStruct = trajStruct(ismember([trajStruct.posture],[1,2,3,4]));
    end
        
%% Choose timestamps 
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    [minNumTimestamps,i] = min(numTimestamps);

    if minNumTimestamps > 10
        minNumTimestamps = 10;
    end
        
%% Update colormaps for dataset
    switch dataset
        case {'E20200316'}
            pcmap = flip(pcmap);
        case {'N20171215'}     
            pcmap = vertcat(pcmap(3,:),pcmap(1,:),pcmap(5,:));
        case {'E20210901'}
            pcmap = vertcat(pcmap(1,:),pcmap(3,:),pcmap(3,:),pcmap(5,:),pcmap(4,:));
            pcmap = vertcat(pcmap(4,:),pcmap(5,:),pcmap(5,:),pcmap(1,:),pcmap(3,:));
        case {'R20201020'}
            pcmap = vertcat(pcmap(2,:),pcmap(1,:));
    end

%% Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);
%     plotTargetList = [1,5];
%     if strcmpi(dataset,'E20200116')
%         plotTargetList = [3,7];
%     end
    
%% Get PCA dimensions
	%Form X, containing trial-averaged data for each condition
    X = NaN(minNumTimestamps,numTargets,numPostures,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            X(:,targetInd,postureInd,:) = traj(1:minNumTimestamps,:); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end

    %Do PCA on all condition averages
    allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);

    % Do marginalizations of X (Xdims: 1=time, 2=target, 3=posture, 4=channel)
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    tSig = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
    pSig = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;

    pSigReshape = reshape(squeeze(pSig),[numPostures*minNumTimestamps,numChannels]);
    tSigReshape = reshape(squeeze(tSig),[numTargets*minNumTimestamps,numChannels]);
    %Compute P and T projections
    [posturePCA,~,~,~,~,~] = pca(pSigReshape); 
    [targetPCA,~,~,~,~,~] = pca(tSigReshape); 

%% Get LDA dimensions
    %Peform LDA by posture on all data 
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = postureList
       obsStruct(structInd).label = posture;
       allObs = [];
       tempTrajStruct = trajStruct([trajStruct.posture]==posture);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
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
    [postureLDA] = doLDA(obsStruct);

    %Peform LDA by target on all data
%     %Use last 3 samples before target acquisition
%     numSamplesToUse = 3;
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for target = targetList
       obsStruct(structInd).label = target;
       allObs = [];
       tempTrajStruct = trajStruct([trajStruct.target]==target);
       timeToUse = 125;
       numOnEitherSide = 1;
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
               traj = tempTrajStruct(i).allSmoothFR(j).traj;
               %[~,timeToUseInd] = min(abs(timestamps-timeToUse));
               timeToUseInd = size(traj,1)-1;
               traj = traj(timeToUseInd-numOnEitherSide:timeToUseInd+numOnEitherSide,:);
               allObs = vertcat(allObs,traj);
%                    if size(traj,1) > minNumTimestamps
%                         traj = traj(1:minNumTimestamps,:);
%                    end
               allObs = vertcat(allObs,traj);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    [targetLDA] = doLDA(obsStruct);
        
%% Orthonormalize; Flip posture axis if necessary 
    %Orthonormalize posture and goal dims
    [PTOrthLDA,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]); PTOrthLDA = PTOrthLDA(:,1:3);
    [PTOrthPCA,~] = qr([posturePCA(:,1),targetPCA(:,1:2)]); PTOrthPCA = PTOrthPCA(:,1:3);

% 
%      %Flip Posture Dim so that posture 1 is on the left
%     p1TrajMean = mean(trajStruct([trajStruct.posture]==1 & [trajStruct.target]==1).avgSmoothFR.traj,1); 
%     p2TrajMean = mean(trajStruct([trajStruct.posture]==max(postureList) & [trajStruct.target]==1).avgSmoothFR.traj,1); 
%     p1TrajCoord = p1TrajMean*postureLDA(:,1);
%     p2TrajCoord = p2TrajMean*postureLDA(:,1);
%     if p1TrajCoord > p2TrajCoord
%         postureLDA(:,1) = -1.*postureLDA(:,1);
%     end
%     p1TrajCoord = p1TrajMean*PTOrth(:,1);
%     p2TrajCoord = p2TrajMean*PTOrth(:,1);
%     if p1TrajCoord > p2TrajCoord
%         PTOrth(:,1) = -1.*PTOrth(:,1);
%     end
% 
%     %Flip C dim if necessary
%     traj = trajStruct(1).avgSmoothFR.traj;
%     CPT = trajStruct(1).avgSmoothFR.traj*CPTOrth;
% %         CP = trajStruct(1).avgSmoothFR.traj*CPOrth;
%     if CPT(end,3) < CPT(1,3)
%         CPTOrth(:,3) = -1.*CPTOrth(:,3);
%     end
% %         if CP(end,3) < CP(1,3)
% %             CPOrth(:,3) = -1.*CPOrth(:,3);
% %         end

%% Add to trajStruct; Compute VAF
    %Add Projections to trajStruct
    totalVar = trace(cov(allTraj));
    for i = 1:size(trajStruct,2)
        %Add average traces to trajStruct
        trajStruct(i).PTOrthLDA.traj = trajStruct(i).avgSmoothFR.traj*PTOrthLDA;
        trajStruct(i).PTOrthPCA.traj = trajStruct(i).avgSmoothFR.traj*PTOrthPCA;
        %Get VAF
        trajStruct(i).PTOrthLDA.VAF =  100.*(diag(cov(allTraj*PTOrthLDA))')./totalVar;
        trajStruct(i).PTOrthPCA.VAF =  100.*(diag(cov(allTraj*PTOrthPCA))')./totalVar;
    end
        
%% Plot 
    %Parameters
    fs = 14;
    timePts = 1:minNumTimestamps;
    
    %Plot - PCA
    figure; hold on      
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrthPCA.traj(timePts,:); 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end   
    grid on
    xticklabels({}); yticklabels({}); zticklabels({}); 
    VAF = round(trajStruct(1).PTOrthPCA.VAF);
    xlabel(['Goal Dim 1 (',num2str(VAF(xDim)),'%)'])
    ylabel(['Goal Dim 2 (',num2str(VAF(yDim)),'%)']) 
    zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
    view([20 10]); set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    if saveFig
        saveas(gcf,fullfile(saveDir,[dataset,'_PTOrth.fig']));
    end
    
    %Plot - LDA
    figure; hold on      
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrthLDA.traj(timePts,:); 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end   
    grid on
    xticklabels({}); yticklabels({}); zticklabels({}); 
    VAF = round(trajStruct(1).PTOrthLDA.VAF);
    xlabel(['Goal Dim 1 (',num2str(VAF(xDim)),'%)'])
    ylabel(['Goal Dim 2 (',num2str(VAF(yDim)),'%)']) 
    zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
    view([20 10]); set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    if saveFig
        saveas(gcf,fullfile(saveDir,[dataset,'_PTOrth.fig']));
    end
    

%% Local functions   
     %Local function for performing LDA - no balancing 
     function [LDAproj] = doLDA(obsStruct)
        numClasses = size(obsStruct,2);
        numDims = size(obsStruct(1).allObs,2);
        obs = []; labels = [];
        for i = 1:numel(obsStruct)
            obs = vertcat(obs,obsStruct(i).allObs);
            labels = vertcat(labels,obsStruct(i).label*ones(size(obsStruct(i).allObs,1),1));
        end
        LDAproj = fisherLDA(obs, labels);
        [LDAproj,~] = qr(LDAproj);
        LDAproj = LDAproj(:,1:numClasses-1);
     end 

%      %LDA - with balancing
% 	function [LDAproj] = doLDA(obsStruct)
%         %Preallocate
%         minNumObs = min([obsStruct.numObs]);
%         numClasses = size(obsStruct,2);
%         numDims = size(obsStruct(1).allObs,2);
%         obs = NaN(minNumObs*numClasses,numDims); labels =  NaN(minNumObs*numClasses,1); 
% 
%         %Fill
%         k = 1;
%         for i = 1:size(obsStruct,2)
%             totalClassObs = size(obsStruct(i).allObs,1);
%             obs(k:k+minNumObs-1,:) = obsStruct(i).allObs(randsample(totalClassObs,minNumObs),:);
%             labels(k:k+minNumObs-1,1) = ones(minNumObs,1).*obsStruct(i).label;
%             k = k+minNumObs;
%         end
% 
%         LDAproj = fisherLDA(obs, labels);
%         [LDAproj,~] = qr(LDAproj);
%         LDAproj = LDAproj(:,1:numClasses-1);
%      end 