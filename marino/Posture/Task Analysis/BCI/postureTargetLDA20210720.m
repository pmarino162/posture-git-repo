clear; clc; clf; close all;

%% Load, Preprocess, and Label Data    
%     date = '20210805';
%     [Data] = loadEarlData20210720;
%     [Data] = loadEarlData20210719;
%     [Data] = loadEarlData20210809;
    [Data] = loadEarlData20210828;
    Data = Data([Data.trialStatus]==1);
    
%% Exlude trials that are too long
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
%     trialInclStates(1).trialName = {'CenterOut_ForceBar_BC'};
%         trialInclStates(1).inclStates = {'BC Freeze','Touch Bar BC'}; 
%         trialInclStates(1).inclOccurrence = {'first','first'};
%     trialInclStates(2).trialName = {'Posture Device Active Movements'};
%         trialInclStates(2).inclStates = {'No Cheating','Target Acquire','Target Hold'};
%         trialInclStates(2).inclOccurrence = {'last','last','last'};
    trialInclStates(2).trialName = {'DelayedCenterOut20210828'};
        trialInclStates(2).inclStates = {'No Cheating','Target Acquire'};
        trialInclStates(2).inclOccurrence = {'last','last'};
    [Data] = excludeLengths(Data,trialInclStates);

%% Run GPFA On States of Interest; Save projection to Spikes Field
    %Info for getStatesTraj    
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Center Target','Step 1'};
        trialInclStates(1).inclOccurrence = {'last','last'};
        trialInclStates(1).addTimeToBeginning = {-100,0};
        trialInclStates(1).addTimeToEnd = {0,100};
%     trialInclStates(1).trialName = {'CenterOut_ForceBar_BC'};
%         trialInclStates(1).inclStates = {'BC Freeze','Touch Bar BC'}; 
%         trialInclStates(1).inclOccurrence = {'first','first'};
%     trialInclStates(2).trialName = {'Posture Device Active Movements'};
    trialInclStates(2).trialName = {'DelayedCenterOut20210828'};
%         trialInclStates(2).inclStates = {'No Cheating','Target Acquire','Target Hold'};
%         trialInclStates(2).inclOccurrence = {'last','last','last'};
%         trialInclStates(2).addTimeToBeginning = {-100,0,0};
%         trialInclStates(2).addTimeToEnd = {0,0,100};
        
        trialInclStates(2).inclStates = {'Center Hold','Delay','No Cheating','Target Acquire','Target Hold'};
        trialInclStates(2).inclOccurrence = {'last','last','last','last','last'};
        trialInclStates(2).addTimeToBeginning = {0,0,0,0,0};
        trialInclStates(2).addTimeToEnd = {0,0,0,0,0};
        

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
    condFields = {{'task','conditionData','taskID'},{'posture','conditionData','postureID'},{'target','targetData','targetID'}};
    trajFields = {'GPFA','marker'};
    
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {0};
        trialInclStates(1).addTimeToEnd = {0};
%     trialInclStates(1).trialName = {'CenterOut_ForceBar_BC'};
%         trialInclStates(1).inclStates = {'BC Freeze','Touch Bar BC'}; 
%         trialInclStates(1).inclOccurrence = {'first','first'};
%         trialInclStates(2).addTimeToBeginning = {0,0};
%         trialInclStates(2).addTimeToEnd = {0,0};  
%     trialInclStates(2).trialName = {'Posture Device Active Movements'};
    trialInclStates(2).trialName = {'DelayedCenterOut20210828'};
%         trialInclStates(2).inclStates = {'No Cheating','Target Acquire','Target Hold'};
%         trialInclStates(2).inclOccurrence = {'last','last','last'};
%         trialInclStates(2).addTimeToBeginning = {0,0,0};
%         trialInclStates(2).addTimeToEnd = {0,0,0};    
        trialInclStates(2).inclStates = {'No Cheating','Target Acquire'};
        trialInclStates(2).inclOccurrence = {'last','last'};
        trialInclStates(2).addTimeToBeginning = {0,0};
        trialInclStates(2).addTimeToEnd = {0,0};    
        
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
            
%% Get number of BC targets and postures
    task = 1;
    tempTrajStruct = trajStruct([trajStruct.task]==task);
    postureList = unique([tempTrajStruct.posture]);
    targetList = unique([tempTrajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    
%% Use BC data to perform posture, target, and PxT LDA
    BCTrajStruct = trajStruct([trajStruct.task]==1);
    allTraj = []; allPostureLabels = []; allTargetLabels = []; allPxTLabels = [];
    for i = 1:size(BCTrajStruct,2)
       target = BCTrajStruct(i).target;
       posture = BCTrajStruct(i).posture;
       numTraj = size(BCTrajStruct(i).allGPFA,2);
       for j = 1:numTraj
            %Trajectory
            traj = BCTrajStruct(i).allGPFA(j).traj;
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
    postureLDA = postureLDA(:,1:numPostures-1);
    
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    targetLDA = targetLDA(:,1:numTargets-1);
%     
    [postTargOrth,~] = qr([postureLDA(:,1:2),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:4);
    
%     [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
%     postTargOrth = postTargOrth(:,1:4);
    
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
    cmap = cmap([1,3,5],:);
%     cmap = cmap([1,2,3,5],:);
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Decomposition Analysis\Earl20200317_Decomposition Analysis\';
    
    for task = 1:2      
        tempTrajStruct = trajStruct([trajStruct.task]==task);
        postureList = unique([tempTrajStruct.posture]);
        targetList = unique([tempTrajStruct.target]);
        
        %Top 3 GPFA - Actual
        figure
        xDim = 1; yDim =2; zDim = 3;
        for posture = postureList
            for target = targetList
                traj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgGPFA.traj;
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
            end
        end
        xlabel(['GPFA ',num2str(xDim)]); ylabel(['GPFA ',num2str(yDim)]); zlabel(['GPFA ',num2str(zDim)])
        grid on
        titleStr = ['Top 3 GPFA']; title(titleStr);
        %saveas(gcf,[dirStr,titleStr,'.fig'])

        %Posture/Target Orth Projection - Actual
        figure
        xDim = 1; yDim = 3; zDim = 4;
        for posture = postureList
            for target = targetList
                postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgPostTargOrth;
                plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
                plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            end
        end
        xlabel(['Posture LDA 1']); ylabel(['Target LDA 1']); zlabel(['Target LDA 2'])
        grid on
        axis equal
        %saveas(gcf,[dirStr,titleStr,'.fig'])     

        figure
        xDim = 2; yDim = 3; zDim = 4;
        for posture = postureList
            for target = targetList
                postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgPostTargOrth;
                plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
                plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
            end
        end
        xlabel(['Posture LDA 2']); ylabel(['Target LDA 1']); zlabel(['Target LDA 2'])
        grid on
        axis equal
        
        figure
        xDim = 1; yDim = 2; zDim = 3;
        for posture = postureList
            for target = targetList
                postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgPostTargOrth;
                plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
                plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
            end
        end
        xlabel(['Posture LDA 1']); ylabel(['Posture LDA 2']); zlabel(['Target LDA 1'])
        grid on
        axis equal
        
        
        %PxT Projection - Actual
        figure
        xDim = 1; yDim = 2; zDim = 3;
        for posture = postureList
            for target = targetList
                pxtProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgPxTLDA;
                plot3(pxtProj(:,xDim),pxtProj(:,yDim),pxtProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
            end
        end
        xlabel(['PxT LDA ',num2str(xDim)]); ylabel(['PxT LDA ',num2str(yDim)]); zlabel(['PxT LDA ',num2str(zDim)])
        grid on
        axis equal
        titleStr = ['PxT']; title(titleStr);
    
    end

    
    %% Posture, Target, and Null vs Time

    for task = 1:2     
        f = figure; f.Position = [10 10 900 300];
        [ha, pos] = tight_subplot(2,5,0.15,0.15,0.05);
        tempTrajStruct = trajStruct([trajStruct.task]==task);
        postureList = unique([tempTrajStruct.posture]);
        targetList = [1,5];
    
        for posture = postureList
            targetInd = 0;
            for target = targetList
                targetInd = targetInd + 1;
                postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgPostTargOrth;
                nullProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgNull;
                time = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgGPFA.timestamps;
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
        sgtitle(['Task ',num2str(task)])
        for sp = 1:10
            axes(ha(sp)); ax = gca; ax.TickDir = 'out';
            xticks([-100 0 250]); xticklabels({'-100',0,'250'})
            ylimits = ylim; ylimits = [round(ylimits(1),1), round(ylimits(2),1)];
%             ylimits = [-2.5, 2.5];
            ax.YLim = ylimits;
            yticks([ylimits(1) ylimits(2)]); yticklabels({num2str(ylimits(1)), num2str(ylimits(2))})
            xlabel('time after go cue (ms)')
            switch sp
                case 1
                    ylabel(['Posture LDA 1 (a.u.)'])
                case 2
                    ylabel(['Posture LDA 2 (a.u.)'])
                case 3
                    ylabel(['Target LDA 1 (a.u.)'])
                case 4
                    ylabel(['Target LDA 2 (a.u.)'])
               
                otherwise
                    ylabel(['Null ',num2str(sp-3),' (a.u.)'])
            end
            set(gca,'fontname','arial')
        end
% 
%     
%         % GPFA vs time
%         f = figure; f.Position = [10 10 1000 400];
%         [ha, pos] = tight_subplot(2,6,0.1,0.1,0.05);
%         for posture = postureList
%             for target = targetList
%                 GPFA = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgGPFA.traj;
%                 time = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgGPFA.timestamps;
%                 for sp = 1:12
%                    axes(ha(sp));
%                    hold on
%                    if target == targetList(1)
%                         plot(time,GPFA(:,sp),'Color',cmap(posture,:),'LineWidth',2); 
%                    elseif target == targetList(2)
%                        plot(time,GPFA(:,sp),':','Color',cmap(posture,:),'LineWidth',2);
%                    end
%                 end
%             end
%         end

    end
    
%% Compare time courses to marker time courses 
%     condFields = {{'task','conditionData','taskID'},{'posture','conditionData','postureID'},{'target','targetData','targetID'}};
%     trajFields = {'GPFA','marker'};
%     
%     conditionData = [Data.conditionData];
%     taskID = [conditionData.taskID];
%     task2Data = Data(taskID==2);
%     
%  
%     trialInclStates(2).trialName = {'Posture Device Active Movements'};
%         trialInclStates(2).inclStates = {'No Cheating','Target Acquire','Target Hold'};
%         trialInclStates(2).inclOccurrence = {'last','last','last'};
%         trialInclStates(2).addTimeToBeginning = {0,0,0};
%         trialInclStates(2).addTimeToEnd = {0,0,0};    
%         
%     markerTrajStruct = getTrajStruct(task2Data,condFields,trajFields,trialInclStates);
% 




    task = 2;     
    f = figure; f.Position = [10 10 500 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.05,0.15);
    tempTrajStruct = trajStruct([trajStruct.task]==task);
    postureList = unique([tempTrajStruct.posture]);
    targetList = [1,5];
    
    for posture = 2
        targetInd = 0;
        for target = targetList
            targetInd = targetInd + 1;
            postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgPostTargOrth;
            nullProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgNull;
            time = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgGPFA.timestamps;
            markerTraj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgMarker.traj;
            markerTime = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgMarker.timestamps;
            for sp = 1:4
               axes(ha(sp));
               hold on
               if targetInd == 1
                    plot(time,postTargOrthProj(:,sp),'Color',cmap(posture,:),'LineWidth',2); 
               elseif targetInd == 2
                   plot(time,postTargOrthProj(:,sp),':','Color',cmap(posture,:),'LineWidth',2);
               end
            end
            for sp = 5:6
               axes(ha(sp));
               hold on
               if targetInd == 1
                    plot(markerTime,markerTraj(:,sp-4),'Color',cmap(posture,:),'LineWidth',2);   
               elseif targetInd == 2
                   plot(markerTime,markerTraj(:,sp-4),':','Color',cmap(posture,:),'LineWidth',2);
               end
            end
        end
    end
    sgtitle(['Task ',num2str(task)])
        
        
    for sp = 1:6
        axes(ha(sp)); ax = gca; ax.TickDir = 'out';
        xticks([-100 0 250]); xticklabels({'-100',0,'250'})
        ylimits = ylim; ylimits = [round(ylimits(1),1), round(ylimits(2),1)];
%             ylimits = [-2.5, 2.5];
        ax.YLim = ylimits;
        yticks([ylimits(1) ylimits(2)]); yticklabels({num2str(ylimits(1)), num2str(ylimits(2))})
        xlabel('time after go cue (ms)')
        switch sp
            case 1
                ylabel(['Posture LDA 1 (a.u.)'])
            case 2
                ylabel(['Posture LDA 2 (a.u.)'])
            case 3
                ylabel(['Target LDA 1 (a.u.)'])
            case 4
                ylabel(['Target LDA 2 (a.u.)'])

            otherwise
                ylabel(['Null ',num2str(sp-3),' (a.u.)'])
        end
        set(gca,'fontname','arial')
    end
% 
