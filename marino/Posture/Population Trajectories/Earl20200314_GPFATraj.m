clear; clc; clf; close all;

%% Load, Preprocess, and Label Data    
    %3/14/2020
    date = '20200314';
    exclCh =  [3 12 16 20 31 73 85 92 94];
    [Data,I30GPFAParams,I30DecoderParams,...
     N00GPFAParams,N00DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200314;
 
%% Run GPFA On States of Interest; Save projection to Spikes Field
    %Info for getStatesTraj    
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
        trialInclStates(1).inclStates = {'Touch Bar Hold','Step 1 Freeze','Step 1'};
        trialInclStates(1).inclOccurrence = {'last','last','last'};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(2).inclStates = {'Touch Bar Hold','Reach'};
        trialInclStates(2).inclOccurrence = {'last','first'};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {'Center Hold','Target'};
        trialInclStates(3).inclOccurrence = {'last','last'};

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
    taskData = [Data.taskData];
    reachingData = Data(strcmpi({taskData.task},'HC'));
    condFields = {{'posture','postureData','postureID'},{'target','targetData','targetID'}};
    trajFields = {'GPFA'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(1).inclStates = {'Touch Bar Hold','Reach'};
        trialInclStates(1).inclOccurrence = {'last','first'};
        trialInclStates(1).addTimeToBeginning = {0};
        trialInclStates(1).addTimeToEnd = {0};
    trajStruct = getTrajStruct(reachingData,condFields,trajFields,trialInclStates);
    
%% Use all data to perform posture, target, and PxT LDA
    allTraj = []; allPostureLabels = []; allTargetLabels = []; allPxTLabels = [];
    for i = 1:size(trajStruct,2)
       posture = trajStruct(i).posture;
       numTraj = size(trajStruct(i).allGPFA,2);
       for j = 1:numTraj
            %Trajectory
            traj = trajStruct(i).allGPFA(j).traj;
            numSteps = size(traj,1);
            allTraj = vertcat(allTraj,traj);
            %Labels
            trajPostureLabel = ones(numSteps,1).*posture;
            allPostureLabels = vertcat(allPostureLabels,trajPostureLabel);
        end
    end
    postureLDA = fisherLDA(allTraj, allPostureLabels);
    [postureLDA,~] = qr(postureLDA);
    postureLDA = postureLDA(:,1:2);
    
    nullSpace = null(postureLDA');
    nullProj = allTraj*nullSpace;
    nullSpacePC = pca(nullProj);
    nullSpace = nullSpace*nullSpacePC;
        
        

%% Add LDA Projections to Traj Struct and Model Struct
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgPostureLDA = trajStruct(i).avgGPFA.traj*postureLDA;
%        trajStruct(i).avgTargetLDA = trajStruct(i).avgGPFA.traj*targetLDA;
%        trajStruct(i).avgPostTargOrth = trajStruct(i).avgGPFA.traj*postTargOrth;
%        trajStruct(i).avgPxTLDA = trajStruct(i).avgGPFA.traj*pxtLDA;
       trajStruct(i).avgNull = trajStruct(i).avgGPFA.traj*nullSpace;
    end

%% Visualize Model and Actual Trajectories 
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Decomposition Analysis\Earl20200317_Decomposition Analysis\';
    
    %Top 3 GPFA - Actual
    figure
    xDim = 1; yDim =2; zDim = 3;
    for posture = 1:5
        for target = [4,5]
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