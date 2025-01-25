clear; clc; clf; close all;

%% Choose Task
%     task = 'BCI';
%     task = 'Iso';
    task = 'Reach';
    
%% Setup save dir   
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20210930 - first model\';

%% Load Data 
    binWidth = 25;
    switch task
        case 'BCI'    
            [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
             N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
             E30GPFAParams,E30DecoderParams] = loadEarlData20200317(binWidth);
        case 'Iso'
            [Data] = loadEarlData20200116(binWidth);
        case 'Reach'
            load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
            exclCh =  [44 87 88 77 78 71 67 69 118];
            getSorts = false;
            Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'binWidth',binWidth,'exclZero',false,'getMarker',true,'getKin',true);
            [Data] = cleanData20210706(Data);
            Data = Data([Data.trialStatus]==1);
            [Data,postureIDs] = labelPostures20210706(Data);
            allTrialPostures = [Data.conditionData];
            allTrialPostures = [allTrialPostures.postureID];
            keepPostures = [1,2,4,5];
            Data = Data(ismember(allTrialPostures,keepPostures));
    end

%% Get Traj Struct   
    trajFields = {'allChannelSmoothedFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    switch task
        case 'BCI'
            condFields = {{'posture','conditionData','postureID'},{'target','targetData','target1ID'}};
            trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
            trialInclStates(1).inclStates = {'Step 1'};
            trialInclStates(1).inclOccurrence = {'last'};
                trialInclStates(1).addTimeToBeginning = {0};
                trialInclStates(1).addTimeToEnd = {0};   
        case 'Iso'
            condFields = {{'posture','conditionData','postureID'},{'target','targetData','targetID'}};
            trialInclStates(1).trialName = {'IsometricForce_1D'};
            trialInclStates(1).inclStates = {'Target'};
            trialInclStates(1).inclOccurrence = {'last'};
                trialInclStates(1).addTimeToBeginning = {0};
                trialInclStates(1).addTimeToEnd = {0}; 
        case 'Reach'
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).trialName = {'GridReaching'};
            trialInclStates(1).inclStates = {'Target Acquire'};
            trialInclStates(1).inclOccurrence = {'first'};
                trialInclStates(1).addTimeToBeginning = {0};
                trialInclStates(1).addTimeToEnd = {0}; 
    end
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,'matchConditions',true);

%% Trim average trajectories to length of shortest one
    avgTrajLengths = nan(1,size(trajStruct,2));
    for i = 1:size(trajStruct,2)
        avgTrajLengths(i) = size(trajStruct(i).avgAllChannelSmoothedFR.traj,1);
    end
    minAvgTrajLength = min(avgTrajLengths);
    for i = 1:size(trajStruct,2)
        trajStruct(i).avgAllChannelSmoothedFR.traj = trajStruct(i).avgAllChannelSmoothedFR.traj(1:minAvgTrajLength,:);
        trajStruct(i).avgAllChannelSmoothedFR.timestamps = trajStruct(i).avgAllChannelSmoothedFR.timestamps(1,1:minAvgTrajLength);
    end
    
%% Get Postures, Targets, and Channels
    postureList = unique([trajStruct.posture]);
    targetList = unique([trajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgAllChannelSmoothedFR.traj,2);
    
%% Do PCA on condition averages
    %Vertically concatenate trial averages
    allAvgs = [];
    for i = 1:size(trajStruct,2)
        traj = trajStruct(i).avgAllChannelSmoothedFR.traj;
        allAvgs = vertcat(allAvgs,traj);
    end
    [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);
       
%% Do LDA on all Data
    postureList = unique([trajStruct.posture]);
    targetList = unique([trajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    allTraj = []; allPostureLabels = []; allTargetLabels = [];
    for i = 1:size(trajStruct,2)
       target = trajStruct(i).target;
       posture = trajStruct(i).posture;
       numTraj = size(trajStruct(i).allAllChannelSmoothedFR,2);
       for j = 1:numTraj
            %Trajectory
            traj = trajStruct(i).allAllChannelSmoothedFR(j).traj;
            numSteps = size(traj,1);
            allTraj = vertcat(allTraj,traj);
            %Labels
            trajTargetLabel = ones(numSteps,1).*target;
            allTargetLabels = vertcat(allTargetLabels,trajTargetLabel);
            trajPostureLabel = ones(numSteps,1).*posture;
            allPostureLabels = vertcat(allPostureLabels,trajPostureLabel);
        end
    end

    postureLDA = fisherLDA(allTraj, allPostureLabels);
    [postureLDA,~] = qr(postureLDA);
    postureLDA = postureLDA(:,1);
    
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    targetLDA = targetLDA(:,1:2);
    
    [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:3);
    
%% Add Projections to trajStruct and modelStruct
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgPCA = (trajStruct(i).avgAllChannelSmoothedFR.traj-mu)*coeff;
       trajStruct(i).avgPostTargOrth = trajStruct(i).avgAllChannelSmoothedFR.traj*postTargOrth;
    end
    
%% Get Variance Explained in condition averages
    %Get total variance in condition averages
        totalVar = sum(var(allAvgs));
%         for i = 1:size(trajStruct,2)
%             totalVar = sum(var(trajStruct(i).avgAllChannelSmoothedFR.traj)) + totalVar;
%         end
    %Get variance explained in each PCA dim
        pcaVarExpl = round((var((allAvgs-mu)*coeff)/totalVar)*100);
    %Get variance explained in each LDA dim
        ldaVarExpl = round((var(allAvgs*postTargOrth)/totalVar)*100);
        
%% Plot in PCA 
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    
    if size(postureList,2) == 2
        cmap = cmap([1,5],:);
    end
    %Top 3 PCA - Actual
    figure
    xDim = 1; yDim =2; zDim = 3;
    postureInd = 1;
    for posture = postureList
        for target = targetList
            traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPCA;
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(traj(5,xDim),traj(5,yDim),traj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(traj(10,xDim),traj(10,yDim),traj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
        end
        postureInd = postureInd + 1;
    end
    xlabel(['PC ',num2str(xDim),' (',num2str(pcaVarExpl(xDim)),'%)']); ylabel(['PC ',num2str(yDim),' (',num2str(pcaVarExpl(yDim)),'%)']); zlabel(['PC ',num2str(zDim),' (',num2str(pcaVarExpl(zDim)),'%)']);
    grid on
    titleStr = ['Actual - Top 3 PC']; title(titleStr);
%     saveas(gcf,[dirStr,task,'ActualTop3PC.fig'])


%% Plot in LDA 
    %Post Targ Orth - Actual
    figure
    xDim = 1; yDim =2; zDim = 3;
    postureInd = 1;
    for posture = postureList
        for target = targetList
            traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPostTargOrth;
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(traj(5,xDim),traj(5,yDim),traj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(traj(10,xDim),traj(10,yDim),traj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
        end
        postureInd = postureInd + 1;
    end
    xlabel(['Posture LDA1 (',num2str(ldaVarExpl(xDim)),'%)']); ylabel(['Target LDA1 (',num2str(ldaVarExpl(yDim)),'%)']); zlabel(['Target LDA2 (',num2str(ldaVarExpl(zDim)),'%)']);
    grid on
    titleStr = ['Actual - LDA']; title(titleStr);
%     saveas(gcf,[dirStr,task,'ActualLDA.fig'])
    
%% PC Timecourses
    CI = zeros(size(trajStruct(1).avgPCA,1),size(trajStruct(1).avgPCA,2)); 
    numTraj = size(trajStruct,2);
    for i = 1:numTraj
       CI = trajStruct(i).avgPCA + CI; 
    end
    CI = CI./numTraj;
    
    figure
    for posture = [1,5]
        for target = [1,5]
            for plotInd = 1:10
                subplot(2,5,plotInd)
                traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPCA;
                traj = traj-CI;
                time = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgAllChannelSmoothedFR.timestamps;
                plot(time,traj(:,plotInd),'Color',cmap(posture,:));
                hold on
            end
        end
    end
    
%% LDA TimeCourses 

    figure
    for posture = [1,5]
        for target = [3,7]
            for plotInd = 1:3
                subplot(1,5,plotInd)
                traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPostTargOrth;
%                 traj = traj-CI;
                time = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgAllChannelSmoothedFR.timestamps;
                plot(time,traj(:,plotInd),'Color',cmap(posture,:));
                hold on
            end
        end
    end