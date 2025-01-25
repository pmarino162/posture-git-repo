clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\Figure 3\Projections';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Get trajStruct
    %Load data
    dataset = 'E20200314';
    [Data,zScoreParams] = loadData(dataset);

    %Get trajStruct
    binWidth = 25; kernelStdDev = 25;
    trialInclStates = struct('trialName','','inclStates',[]);
    trajFields = {'zSmoothFR','markerVel','marker'};
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};

    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
        trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',0},{'state','Step 2','first',0}};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

    taskIDs = struct('ID',[],'task','');
    taskIDs(1).ID = 1; taskIDs(1).task = 'BC';
    taskIDs(2).ID = 2; taskIDs(2).task = 'HC';
    taskIDs(3).ID = 3; taskIDs(3).task = 'Iso';
    taskIDs(4).ID = 4; taskIDs(4).task = 'All';
        
%% For each task, get posture signal and space
    resultStruct = struct('task',[],'pSig',[],'pSpace',[]);
    for task = 1:3
        if task == 4
            taskTrajStruct = trajStruct;
        else
            taskTrajStruct = trajStruct([trajStruct.task]==task);
        end
        
        % Get minimum number of timestamps in condition averages
        numTimestamps = [];
        for i = 1:size(taskTrajStruct,2)
            numTimestamps(i) = length(taskTrajStruct(i).avgSmoothFR.timestamps);
        end
        [minNumTimestamps,i] = min(numTimestamps);
               
        if minNumTimestamps > 10
            numPts = 10;
        else
            numPts = minNumTimestamps;
        end
        
        %Get numPostures and numTargets
        postureList = unique([taskTrajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([taskTrajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(taskTrajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(taskTrajStruct,2);
        
        
        [pSig,tSig] = getPandTsig(taskTrajStruct,numPts);
        pSigReshape = reshape(squeeze(pSig),[numPostures*numPts,numChannels]);
        tSigReshape = reshape(squeeze(tSig),[numTargets*numPts,numChannels]);
        [pSpace,~,~,~,explainedP,pSigMu] = pca(pSigReshape); 
        [tSpace,~,~,~,explainedT,tSigMu] = pca(tSigReshape); 
        
        resultStruct(task).task = task;
        resultStruct(task).pSig = pSig;
        resultStruct(task).pSpace = pSpace;
        
    end
    
%% For each task, plot posture signal in each posture space
    figure; hold on;
    plotInd = 1;
    for spaceTask = 1:3
        pSpace = resultStruct([resultStruct.task]==spaceTask).pSpace;
        for sigTask = 1:3
            pSig = resultStruct([resultStruct.task]==sigTask).pSig;
            for posture = 1:3
                proj = squeeze(pSig(:,1,posture,:))*pSpace;
                subplot(3,3,plotInd); hold on
                plot(proj(:,1),proj(:,2),'LineWidth',2,'Color',pcmap(posture,:))
            end
            plotInd = plotInd + 1;
        end        
    end
    
    
%% Local functions
    function [pSig,tSig] = getPandTsig(trajStruct,numPts)
        %Get posture & target lists
        postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
        numManPCs = size(trajStruct(1).avgSmoothFR.traj,2); 

        %Form X
        X = NaN(numPts,numTargets,numPostures,numManPCs);
        postureInd = 1;
        for posture = postureList
            targetInd = 1;
            for target = targetList
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
                X(:,targetInd,postureInd,:) = traj(1:numPts,:); 
                targetInd = targetInd + 1;
            end
            postureInd = postureInd+1;
        end

        % Do marginalizations of X (Xdims: 1=time, 2=target, 3=posture, 4=channel)
        %Condition-invariant
        CIMargOffset = mean(X,[1 2 3],'omitnan');
        CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
        %Posture and Target Traj
        tSig = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
        pSig = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
    end
