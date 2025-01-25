clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20221123\Figures\Supplement\Tuning Curves';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = vertcat(orli,customRainbow(4:5,:));
    tcmap = customRainbow;

%% Run Loop    
    datasetList = {'R20201020'};
       clf; close all
        %{'E20200317','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}
        % Load Data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR','bciCursorPos','bciCursorVel'};
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
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',250}};
            case 'E20211007'
                trialInclStates(1).trialName = {'IsometricForce_1D'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
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
    
%% Get PCs
    % Get minimum number of timestamps in condition averages
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    [minNumTimestamps,i] = min(numTimestamps);

    %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numCh = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);
        
    %PC timecourses        
    %Form X, containing trial-averaged data for each condition
    X = NaN(minNumTimestamps,numTargets,numPostures,numCh);
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
  
    %Do PCA; mean center; project down to manifold
    allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numCh]);
    [allPCs,~,~,~,explained,allMu] = pca(allTraj);
            
%% Plot all cursor trajectories by target
    for posture = 1:2
        figure; hold on
        for target = 1:8
            allPos = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture).allBciCursorPos;
            for i = 1:size(allPos,2)
                traj = allPos(i).traj;
                plot(traj(:,1),traj(:,2),'Color',tcmap(target,:))
            end
        end
    end

%% Compare cursor traj and FR to mean PSTH
    channels = [2,3,4];
    for posture = 1:2
        figure; 
        for target = 1:8
            pos = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture).avgBciCursorPos.traj;
            posTime = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture).avgBciCursorPos.timestamps;
            vel = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture).avgBciCursorVel.traj;
            velTime = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture).avgBciCursorVel.timestamps;
            neural = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture).avgSmoothFR.traj;
            neuralTime = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture).avgSmoothFR.timestamps;
            
            popFR = mean(neural,2);
            
            subplot(8,1,1)
                plot(posTime,pos(:,1),'Color',tcmap(target,:))
                ylabel('X'); hold on
            subplot(8,1,2)
                plot(posTime,pos(:,2),'Color',tcmap(target,:))
                ylabel('Y'); hold on
            subplot(8,1,3)
                plot(velTime,vel(:,1),'Color',tcmap(target,:))
                ylabel('Vx'); hold on
            subplot(8,1,4)
                plot(velTime,vel(:,2),'Color',tcmap(target,:))
                ylabel('Vy'); hold on
            subplot(8,1,5)
                plot(neuralTime,popFR,'Color',tcmap(target,:))
                ylabel('PopFR'); hold on
            subplot(8,1,6)
                plot(neuralTime,neural(:,7),'Color',tcmap(target,:))
                ylabel('Ch'); hold on
            subplot(8,1,7)
                plot(neuralTime,neural(:,8),'Color',tcmap(target,:))
                ylabel('Ch'); hold on
            subplot(8,1,8)
                plot(neuralTime,neural(:,9),'Color',tcmap(target,:))
                ylabel('Ch'); hold on

        end
    end
    
%% For particular trial, plot cursor motion and population FR
trial = 2;
[traj,timestamps] = getStatesTraj20220419(Data(trial),trialInclStates,'zSmoothFR',binWidth,kernelStdDev,'zScoreParams',zScoreParams,'rmSort',[]);
popFR = mean(traj,2)';

pos = Data(trial).Decoder.position;
posTime = Data(trial).Decoder.posTime;
vel = Data(trial).Decoder.velocity;
velTime = Data(trial).Decoder.velTime;
stateTransitions = Data(trial).stateData.stateTransitions;

        reactTime = stateTransitions(2,find(stateTransitions(1,:)==4));
        moveTime = stateTransitions(2,find(stateTransitions(1,:)==5));
        holdTime = stateTransitions(2,find(stateTransitions(1,:)==6));
       
        %Plot decoded pos/vel with state transtions
        figure
            subplot(5,1,1); hold on
                plot(posTime,pos(1,:))
                ax = gca;
                line([reactTime,reactTime],[ax.YLim])
                line([moveTime,moveTime],[ax.YLim])
                line([holdTime,holdTime],[ax.YLim])
                ylabel('x')
            subplot(5,1,2); hold on
                plot(posTime,pos(2,:))
                                ax = gca;
                line([reactTime,reactTime],[ax.YLim])
                line([moveTime,moveTime],[ax.YLim])
                line([holdTime,holdTime],[ax.YLim])
                ylabel('y')
            subplot(5,1,3); hold on
                plot(velTime,vel(1,:))
                                ax = gca;
                line([reactTime,reactTime],[ax.YLim])
                line([moveTime,moveTime],[ax.YLim])
                line([holdTime,holdTime],[ax.YLim])
                ylabel('Vx')
            subplot(5,1,4); hold on
                plot(velTime,vel(2,:))
                                ax = gca;
                line([reactTime,reactTime],[ax.YLim])
                line([moveTime,moveTime],[ax.YLim])
                line([holdTime,holdTime],[ax.YLim])
                ylabel('Vy')
                
            subplot(5,1,5); hold on
                plot(timestamps + reactTime,popFR)
                
                ax = gca;
                ax.XLim = [velTime(1) velTime(end)]
                                line([reactTime,reactTime],[ax.YLim])
                line([moveTime,moveTime],[ax.YLim])
                line([holdTime,holdTime],[ax.YLim])