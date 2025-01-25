clear; clc; clf; close all;

%% Load Data 
    %Earl - 3/17/2020
    date = '20200317';
    exclCh =  [3 12 16 20 31 73 85 92 94];
    binWidth = 1;
    [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
     N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200317(binWidth);
    
%% Run GPFA On States of Interest; Save projection to Spikes Field
    %Info for getStatesTraj
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
        trialInclStates(1).inclStates = {'Step 1'}; trialInclStates(1).inclOccurrence = {'last'};
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
    condFields = {{'posture','postureData','postureID'},{'target','targetData','target1ID'}};
    trajFields = {'GPFA'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
        trialInclStates(1).inclStates = {'Step 1'}; trialInclStates(1).inclOccurrence = {'last'};
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
    
    [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:3);
    
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
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Paper Figures\Trajectories\';

    f = figure; f.Position = [10 10 800 300];
    f.Renderer = 'Painters';
    [ha, pos] = tight_subplot(1,3,0,0.15,0.07);
    for sp = 1:3
        xDim = 1; yDim = 2; zDim = 3;
        axes(ha(sp));
        for posture = 1:5
            for target = 1:8
                postTargOrthProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPostTargOrth;
                plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
            end
        end
        xlabel(['Posture LDA 1']); ylabel(['Target LDA 1']); zlabel(['Target LDA 2'])
        grid on
        axis equal   
        switch sp
            case 1
                
            case 2
                view(90,0)
            case 3
                view(0,90)
        end
        ax = gca;
        ax.TickDir = 'out';
        xticklabelcell = ax.XTickLabel;
        xticklabelcell(2:end-1,:) = {''};
        ax.XTickLabel = xticklabelcell;

        yticklabelcell = ax.YTickLabel;
        yticklabelcell(2:end-1,:) = {''};
        ax.YTickLabel = yticklabelcell;
        
        zticklabelcell = ax.ZTickLabel;
        zticklabelcell(2:end-1,:) = {''};
        ax.ZTickLabel = zticklabelcell;

        set(gca,'fontname','arial')
    
    end
    titleStr = 'Earl_3DTraj';
    saveas(gcf,[dirStr,titleStr,'.svg'])  

    
    %% Posture, Target, and Null vs Time
    close all
    f = figure; f.Position = [10 10 700 300];
    f.Renderer = 'Painters';
    [ha, pos] = tight_subplot(2,3,0.15,0.15,0.05);
    
    
    for posture = [1,5]
        for target = [1,5]
            postTargOrthProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPostTargOrth;
            nullProj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgNull;
            time = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgGPFA.timestamps;
            for sp = 1:3
                axes(ha(sp));
               hold on
               if target == 1
                    plot(time,postTargOrthProj(:,sp),'Color',cmap(posture,:),'LineWidth',2); 
               elseif target == 5
                   plot(time,postTargOrthProj(:,sp),':','Color',cmap(posture,:),'LineWidth',2);
               end
            end
            for sp = 4:5
                axes(ha(sp));
               hold on
               if target == 1
                    plot(time,nullProj(:,sp-3),'Color',cmap(posture,:),'LineWidth',2);   
               elseif target == 5
                   plot(time,nullProj(:,sp-3),':','Color',cmap(posture,:),'LineWidth',2);
               end
            end
        end
    end
    
    axes(ha(1))
    ax = gca;
    ax.TickDir = 'out';
    xticks([-100 0 250])
    xticklabels({'-100',0,'250'})
    ylimits = ylim;
    ylimits = [round(ylimits(1),1), round(ylimits(2),1)];
    yticks([ylimits(1) ylimits(2)])
    yticklabels({num2str(ylimits(1)), num2str(ylimits(2))})
    xlabel('time after go cue (ms)')
    ylabel('Posture LDA 1 (a.u.)')
    set(gca,'fontname','arial')
    
    
    axes(ha(2))
    ax = gca;
    ax.TickDir = 'out';
    xticks([-100 0 250])
    xticklabels({'-100',0,'250'})
    ylimits = ylim;
    ylimits = [round(ylimits(1),1), round(ylimits(2),1)];
    yticks([ylimits(1) ylimits(2)])
    yticklabels({num2str(ylimits(1)), num2str(ylimits(2))})
    xlabel('time after go cue (ms)')
    ylabel('Target LDA 1 (a.u.)')
    set(gca,'fontname','arial')
    
    axes(ha(3))
    ax = gca;
    ax.TickDir = 'out';
    xticks([-100 0 250])
    xticklabels({'-100',0,'250'})
    ylimits = ylim;
    ylimits = [round(ylimits(1),1), round(ylimits(2),1)];
    yticks([ylimits(1) ylimits(2)])
    yticklabels({num2str(ylimits(1)), num2str(ylimits(2))})
    xlabel('time after go cue (ms)')
    ylabel('Target LDA 2 (a.u.)')
    set(gca,'fontname','arial')
    
    axes(ha(4))
    ax = gca;
    ax.TickDir = 'out';
    xticks([-100 0 250])
    xticklabels({'-100',0,'250'})
    ylimits = ylim;
    ylimits = [round(ylimits(1),1), round(ylimits(2),1)];
    yticks([ylimits(1) ylimits(2)])
    yticklabels({num2str(ylimits(1)), num2str(ylimits(2))})
    xlabel('time after go cue (ms)')
    ylabel('Null 1 (a.u.)')
    set(gca,'fontname','arial')
    
    axes(ha(5))
    ax = gca;
    ax.TickDir = 'out';
    xticks([-100 0 250])
    xticklabels({'-100',0,'250'})
    ylimits = ylim;
    ylimits = [round(ylimits(1),1), round(ylimits(2),1)];
    yticks([ylimits(1) ylimits(2)])
    yticklabels({num2str(ylimits(1)), num2str(ylimits(2))})
    xlabel('time after go cue (ms)')
    ylabel('Null 2 (a.u.)')
    set(gca,'fontname','arial')
    
    titleStr = 'Earl_5DTraj';
%     saveas(gcf,[dirStr,titleStr,'.svg'])  

    