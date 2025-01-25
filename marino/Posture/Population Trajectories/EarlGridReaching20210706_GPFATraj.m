clear; clc; clf; close all;

%% Load Data 
    %Earl - 7/06/2021
    date = '20210706';
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
    
%% Clean Data
    [Data] = cleanData20210706(Data);
    
%% Get Data Struct, Keep only successful trials
    exclCh =  [44 87 88 77 78 71 67 69 118];
    exclCh = [];
    getSorts = false;
    Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'exclZero',false,'getMarker',true,'getKin',true);
    Data = Data([Data.trialStatus]==1);
    
%% Label Postures
    %Get list of targets for grid reaching trials
    gridReachData = Data(cellfun(@(x) strcmpi(x,'GridReaching'),{Data.trialName}));
    targetData = [gridReachData.targetData];
    workspaceCenter = cell2mat({targetData.workspaceCenter}');
    uniWorkspaceCenter = unique(workspaceCenter,'rows');
    sortedYs = sort(unique(uniWorkspaceCenter(:,2)),'descend')';
    sortedXs = sort(unique(uniWorkspaceCenter(:,1)))';
    numCol = size(sortedXs,2);
    numTrials = size(gridReachData,2);
    for trial = 1:numTrials
       workspaceCenter = gridReachData(trial).targetData.workspaceCenter(1:2);
       row = find(sortedYs==workspaceCenter(2));
       col = find(sortedXs==workspaceCenter(1));
       gridReachData(trial).conditionData.postureID = (row-1)*numCol + col;
    end
    %Create postureIDs list
    postureIDs = struct('postureID',[],'workspaceCenter',[]);
    for i = 1:size(uniWorkspaceCenter,1)
       workspaceCenter = uniWorkspaceCenter(i,:);
       row = find(sortedYs==workspaceCenter(2));
       col = find(sortedXs==workspaceCenter(1));
       postureID = (row-1)*numCol + col;
       postureIDs(i).postureID = postureID;
       postureIDs(i).workspaceCenter = workspaceCenter;      
    end
    [~,sortInd] = sort([postureIDs.postureID]);
    postureIDs = postureIDs(sortInd);
    
 %% Keep only gridReachingData; 
    Data = gridReachData;
    clearvars gridReachData 
    
%% Run GPFA On States of Interest; Save projection to Spikes Field
    %Info for getStatesTraj
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'GridReaching'};
        trialInclStates(1).inclStates = {'Target Acquire'}; trialInclStates(1).inclOccurrence = {'first'};
        trialInclStates(1).addTimeToBeginning = {-250};
        trialInclStates(1).addTimeToEnd = {250};
    %Info for GPFA
    numDims = 12;
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
    condFields = {{'posture','conditionData','postureID'},{'target','targetData','targetID'}};
    trajFields = {'GPFA'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'GridReaching'};
        trialInclStates(1).inclStates = {'Target Acquire'}; trialInclStates(1).inclOccurrence = {'first'};
        trialInclStates(1).addTimeToBeginning = {-100};
        trialInclStates(1).addTimeToEnd = {0};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
    
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
    
    [postTargOrth,~] = qr([postureLDA(:,1:2),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:4);
    
    nullSpace = null(postTargOrth');
    nullProj = allTraj*nullSpace;
    nullSpacePC = pca(nullProj);
    nullSpace = nullSpace*nullSpacePC;
        
    pxtLDA = fisherLDA(allTraj, allPxTLabels);
    [pxtLDA,~] = qr(pxtLDA);
    pxtLDA = pxtLDA(:,1:numDims-1);
        

%% Add LDA Projections to Traj Struct and Model Struct
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgPostureLDA = trajStruct(i).avgGPFA.traj*postureLDA;
       trajStruct(i).avgTargetLDA = trajStruct(i).avgGPFA.traj*targetLDA;
       trajStruct(i).avgPostTargOrth = trajStruct(i).avgGPFA.traj*postTargOrth;
       trajStruct(i).avgPxTLDA = trajStruct(i).avgGPFA.traj*pxtLDA;
       trajStruct(i).avgNull = trajStruct(i).avgGPFA.traj*nullSpace;
    end

%% Visualize Model and Actual Trajectories 
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Decomposition Analysis\Earl20200317_Decomposition Analysis\';
    
    %Top 3 GPFA - Actual
    figure
    xDim = 4; yDim =5; zDim = 6;
    for posture = [1,5]
        for target = [1,5]
            traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgGPFA.traj;
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['GPFA ',num2str(xDim)]); ylabel(['GPFA ',num2str(yDim)]); zlabel(['GPFA ',num2str(zDim)])
    grid on
    titleStr = ['Actual - Top 3 GPFA']; title(titleStr);
%     saveas(gcf,[dirStr,titleStr,'.fig'])

    %Posture/Target Orth Projection - Actual
    figure
    xDim = 2; yDim = 3; zDim = 1;
    for posture = 1:5
        for target = [4,5]
            postTargOrthProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPostTargOrth;
            plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['Target LDA 1']); ylabel(['Target LDA 2']); zlabel(['Posture LDA 1'])
    grid on
    axis equal
% %     saveas(gcf,[dirStr,titleStr,'.fig'])     

    %PxT Projection - Actual
    figure
    xDim = 1; yDim = 2; zDim = 3;
    for posture = 1:5
        for target = [4,5]
            pxtProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPxTLDA;
            plot3(pxtProj(:,xDim),pxtProj(:,yDim),pxtProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    xlabel(['PxT LDA ',num2str(xDim)]); ylabel(['PxT LDA ',num2str(yDim)]); zlabel(['PxT LDA ',num2str(zDim)])
    grid on
    axis equal
    titleStr = ['Actual - PxT']; title(titleStr);
    


    
    %% Posture, Target, and Null vs Time
    f = figure; f.Position = [10 10 900 300];
    [ha, pos] = tight_subplot(2,5,0.15,0.15,0.05);
    
    for posture = [1,5]
        targetInd = 0;
        for target = [4,6]
            targetInd = targetInd + 1;
            postTargOrthProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPostTargOrth;
            nullProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgNull;
            time = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgGPFA.timestamps;
            for sp = 1:4
               axes(ha(sp));
               hold on
               if targetInd == 1
                    plot(time,postTargOrthProj(:,sp),'Color',cmap(posture,:),'LineWidth',2); 
               elseif targetInd == 2
                   plot(time,postTargOrthProj(:,sp),':','Color',cmap(posture,:),'LineWidth',2);
               end
            end
            for sp = 5:10
                axes(ha(sp));
               hold on
               if targetInd == 1
                    plot(time,nullProj(:,sp-4),'Color',cmap(posture,:),'LineWidth',2);   
               elseif targetInd == 2
                   plot(time,nullProj(:,sp-4),':','Color',cmap(posture,:),'LineWidth',2);
               end
            end
        end
    end
    sgtitle(['Targets ',num2str(target)])
    for sp = 1:10
        axes(ha(sp)); ax = gca; ax.TickDir = 'out';
        xticks([-100 0 250]); xticklabels({'-100',0,'250'})
%         ylimits = ylim; ylimits = [round(ylimits(1),1), round(ylimits(2),1)];
        ylimits = [-2.5, 2.5];
        ax.YLim = ylimits;
        yticks([ylimits(1) ylimits(2)]); yticklabels({num2str(ylimits(1)), num2str(ylimits(2))})
        xlabel('time after go cue (ms)')
        switch sp
            case 1
                ylabel(['Posture LDA 1 (a.u.)'])
            case 2
                ylabel(['Target LDA 1 (a.u.)'])
            case 3
                ylabel(['Target LDA 2 (a.u.)'])
            otherwise
                ylabel(['Null ',num2str(sp-3),' (a.u.)'])
        end
        set(gca,'fontname','arial')
    end

    
%% GPFA vs time
    f = figure; f.Position = [10 10 1000 400];
    [ha, pos] = tight_subplot(2,6,0.1,0.1,0.05);
    for posture = [1,5]
        for target = [4,8]
            GPFA = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgGPFA.traj;
            time = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgGPFA.timestamps;
            for sp = 1:12
               axes(ha(sp));
               hold on
               if target == 4
                    plot(time,GPFA(:,sp),'Color',cmap(posture,:),'LineWidth',2); 
               elseif target == 8
                   plot(time,GPFA(:,sp),':','Color',cmap(posture,:),'LineWidth',2);
               end
            end
        end
    end
    