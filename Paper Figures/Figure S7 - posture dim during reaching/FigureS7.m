clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S7 - posture dim during reaching';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Parameters
    numPCsToKeep = 10;
    regularize = true;
    
%% Run loop for each dataset  
reachDatasetList = {'N20190226','R20200221','E20210706'};

for datasetList = {'E20210706'};%reachDatasetList%{'E20210706'}
    %Load data
    dataset = datasetList{1,1};
    [Data,zScoreParams] = loadData(dataset);
    
    %Get trajStruct
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    trajFields = {'zSmoothFR','marker'};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);                     
    %Get moveTrajStruct - align to move onset
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',300}};
    moveTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);
    %Get holdTrajStruct - align to trial end
    trialInclStates(1).inclStates = {{'state','Success with Reward','first',-200},{'state','Success with Reward','first',0}};
    holdTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);
    % Keep postures 2&7 only for earl reaching
    if strcmpi(dataset,'E20210706')
        trajStruct = trajStruct(ismember([trajStruct.posture],[2,7]));
        moveTrajStruct = moveTrajStruct(ismember([moveTrajStruct.posture],[2,7]));
        holdTrajStruct = holdTrajStruct(ismember([holdTrajStruct.posture],[2,7]));
    end
    %Get trajStructDimensions
    [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);

    %% Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(numPostures);    
    
    %% PCA version
    %Project all data down to top PCs
    [allTraj,Mu] = collectAllAvgTraj(trajStruct);
    [trajStruct,PCs] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep);
    % Fit P and T spaces, components for each group
    [minNumTimestamps] = getMinNumTimestamps(trajStruct);    
    [pSig,tSig] = getPandTsig(trajStruct,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
    pSigReshape = reshape(squeeze(pSig),[numPostures*minNumTimestamps,numPCsToKeep]);
    [pDims,~,~,~,explainedP,pSigMu] = pca(pSigReshape); 
    
    %% dPCA Version
    [W,V,postureDims,targetDims,explVar,additionalVarExpl] = getPostureDPCADims(trajStruct,regularize,postureList,numPostures,targetList,numTargets,numChannels,minNumTimestamps);
    pDPCA = W(:,postureDims);
    pDPCA = -1*pDPCA;
    
    %% Get marker axis 
    targetID = [Data.targetData];
    targetID = [targetID.targetID];
    postureID = [Data.conditionData];
    postureID = [postureID.postureID];
    switch dataset
        case {'N20190226'}
            target4Pos1Loc = [-27.4280, 67.7280];
            target8Pos1Loc = [57.4240, -17.1240];
            target8Pos2Loc = [159.8240, -86.5240];
            targetAxis = (target4Pos1Loc-target8Pos1Loc)./norm(target4Pos1Loc-target8Pos1Loc);
            targetAxisCenterLoc = target8Pos1Loc;
            
        case {'R20200221'}
            target4Pos1Loc = [45.3720, 87.3280];
            target8Pos1Loc = [130.2240, 2.4760];
            target8Pos2Loc = [232.6240, -66.9240];
            targetAxis = (target4Pos1Loc-target8Pos1Loc)./norm(target4Pos1Loc-target8Pos1Loc);
            targetAxisCenterLoc = target8Pos1Loc;

            
        case{'E20210706'}
            target4Pos2Loc = Data(find(targetID==4 & postureID==2,1,'first')).targetData.targetLoc(1,1:2) + ...
                Data(find(targetID==4 & postureID==2,1,'first')).targetData.workspaceCenter(1,1:2);
            target8Pos7Loc = Data(find(targetID==8 & postureID==7,1,'first')).targetData.targetLoc(1,1:2) + ...
                Data(find(targetID==8 & postureID==7,1,'first')).targetData.workspaceCenter(1,1:2);
            targetAxis = (target4Pos2Loc-target8Pos7Loc)./norm(target4Pos2Loc-target8Pos7Loc);
            targetAxisCenterLoc = target8Pos7Loc + (target4Pos2Loc-target8Pos7Loc)/2;
    end
    
%% Project data onto relevant axes
    for i = 1:size(trajStruct,2)
        %Add average traces to trajStruct
        trajStruct(i).pPCA.traj = trajStruct(i).avgPCA.traj*pDims;       
        
        moveTrajStruct(i).avgPCA.traj = (moveTrajStruct(i).avgZSmoothFR.traj-Mu)*PCs;
        moveTrajStruct(i).pPCA.traj = moveTrajStruct(i).avgPCA.traj(:,1:numPCsToKeep)*pDims;
        moveTrajStruct(i).pPCA.timestamps = moveTrajStruct(i).avgZSmoothFR.timestamps;    
        moveTrajStruct(i).pDPCA.traj = moveTrajStruct(i).avgZSmoothFR.traj*pDPCA;
        moveTrajStruct(i).pDPCA.timestamps = moveTrajStruct(i).avgZSmoothFR.timestamps;      
        moveTrajStruct(i).targetAxis.traj = (moveTrajStruct(i).avgMarker.traj(:,1:2)-targetAxisCenterLoc(:,1:2))*targetAxis';
        moveTrajStruct(i).targetAxis.timestamps = moveTrajStruct(i).avgMarker.timestamps;
             
        holdTrajStruct(i).avgPCA.traj = (holdTrajStruct(i).avgZSmoothFR.traj-Mu)*PCs;
        holdTrajStruct(i).pPCA.traj = holdTrajStruct(i).avgPCA.traj(:,1:numPCsToKeep)*pDims;
        holdTrajStruct(i).pPCA.timestamps = holdTrajStruct(i).avgZSmoothFR.timestamps;     
        holdTrajStruct(i).pDPCA.traj = holdTrajStruct(i).avgZSmoothFR.traj*pDPCA;
        holdTrajStruct(i).pDPCA.timestamps = holdTrajStruct(i).avgZSmoothFR.timestamps;      
        holdTrajStruct(i).targetAxis.traj = (holdTrajStruct(i).avgMarker.traj(:,1:2)-targetAxisCenterLoc(:,1:2))*targetAxis';
        holdTrajStruct(i).targetAxis.timestamps = holdTrajStruct(i).avgMarker.timestamps;
        
        %Add all traces to struct
        for j = 1:size(trajStruct(i).allZSmoothFR,2)
            moveTrajStruct(i).allPCA(j).traj = (moveTrajStruct(i).allZSmoothFR(j).traj-Mu)*PCs;
            moveTrajStruct(i).allPPCA(j).traj = moveTrajStruct(i).allPCA(j).traj(:,1:numPCsToKeep)*pDims;
            moveTrajStruct(i).allPPCA(j).timestamps = moveTrajStruct(i).allZSmoothFR(j).timestamps;
            moveTrajStruct(i).allTargetAxis(j).traj = (moveTrajStruct(i).allMarker(j).traj(:,1:2)-targetAxisCenterLoc(:,1:2))*targetAxis';
            moveTrajStruct(i).allTargetAxis(j).timestamps = moveTrajStruct(i).allMarker(j).timestamps;         
            moveTrajStruct(i).allDPCA(j).traj = moveTrajStruct(i).allZSmoothFR(j).traj*pDPCA;
            moveTrajStruct(i).allDPCA(j).timestamps = moveTrajStruct(i).allZSmoothFR(j).timestamps;
        end
    end

%% Plot

switch dataset
    case {'N20190226','R20200221'}
        plotPostureList = [1,2];
        plotTargetList = [4,8];
    case {'E20210706'}
        plotPostureList = [2,7];
        plotTargetList = [4,8];
end

traceFigScale = .3;
traceFs = 6;
traceLineWidth = 1;
corrFigScale = .2;
corrFs = 6;
corrDotSize = 1;

timeOffset = 550;
%Posture axis projection
f = figure; hold on
f.Position = [200 200 f.Position(3)*traceFigScale f.Position(4)*traceFigScale*1.5]; 

postureInd = 1;
for posture = plotPostureList
    for target = plotTargetList
        %Posture Axis (Neural Activity)
        subplot(2,1,1);hold on
        moveTraj = moveTrajStruct([moveTrajStruct.posture]==posture & [moveTrajStruct.target]==target).pPCA.traj(:,1);
        moveTime = moveTrajStruct([moveTrajStruct.posture]==posture & [moveTrajStruct.target]==target).pPCA.timestamps;        
        holdTraj = holdTrajStruct([holdTrajStruct.posture]==posture & [holdTrajStruct.target]==target).pPCA.traj(:,1);
        holdTime = holdTrajStruct([holdTrajStruct.posture]==posture & [holdTrajStruct.target]==target).pPCA.timestamps;
        holdTime = holdTime + timeOffset;        
        if target == plotTargetList(1)
            plot(moveTime,moveTraj,'--','LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
            plot(holdTime,holdTraj,'--','LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
        elseif target == plotTargetList(2)
            plot(moveTime,moveTraj,'LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
            plot(holdTime,holdTraj,'LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
        end          
        
        %Hand (Along target axis)
        subplot(2,1,2);hold on
        moveTraj = moveTrajStruct([moveTrajStruct.posture]==posture & [moveTrajStruct.target]==target).targetAxis.traj(:,1);
        moveTime = moveTrajStruct([moveTrajStruct.posture]==posture & [moveTrajStruct.target]==target).targetAxis.timestamps;
        holdTraj = holdTrajStruct([holdTrajStruct.posture]==posture & [holdTrajStruct.target]==target).targetAxis.traj(:,1);
        holdTime = holdTrajStruct([holdTrajStruct.posture]==posture & [holdTrajStruct.target]==target).targetAxis.timestamps;
        holdTime = holdTime + timeOffset;      
        if target == plotTargetList(1)
            plot(moveTime,moveTraj,'--','LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
            plot(holdTime,holdTraj,'--','LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
        elseif target == plotTargetList(2)
            plot(moveTime,moveTraj,'LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
            plot(holdTime,holdTraj,'LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
        end        
        
    end
    postureInd = postureInd + 1;
end

subplot(2,1,1)
    xticks([-200, 0, 300, 350, 550])
    xticklabels({})
    ylabel(['Posture',newline,'Dim 1 (a.u.)'])
    ax = gca;
    rectangle('Position',[0,ax.YLim(1),300,ax.YLim(2)-ax.YLim(1)], 'FaceColor',[0.7 0.7 0.7 .25], 'EdgeColor',[1 1 1 0])
    ax.YTick = [ax.YTick(1) ax.YTick(end)];
    ax.TickDir = 'out';
    ax.FontName = 'arial';
    ax.FontSize = traceFs;
    
subplot(2,1,2)
    xticks([-200, 0, 300, 350, 550])
    xticklabels({'Move-200','Move','Move+300','Trial End-200','Trial End'})
    xtickangle(45) 
    xlabel('time (ms)')
    ylabel(['Hand',newline,'Pos (mm)'])
    ax = gca;
    ax.YTick = [ax.YTick(1) ax.YTick(end)];
    ylim = ax.YLim;
    rectangle('Position',[0,ax.YLim(1),300,ax.YLim(2)-ax.YLim(1)], 'FaceColor',[0.7 0.7 0.7 .25], 'EdgeColor',[1 1 1 0])
    ax.YLim = ylim;
    ax.TickDir = 'out';
    ax.TickDir = 'out';
    ax.FontName = 'arial';
    ax.FontSize = traceFs;
    
    if saveFig
        saveas(gcf,fullfile(saveDir,[dataset,'_PCA_Traces.svg']));
    end
    
    
%Posture axis projection
f = figure; hold on
f.Position = [200 200 f.Position(3)*traceFigScale f.Position(4)*traceFigScale*1.5]; 
postureInd = 1;
for posture = plotPostureList
    for target = plotTargetList
        %Posture Axis (Neural Activity)
        subplot(2,1,1);hold on
        moveTraj = moveTrajStruct([moveTrajStruct.posture]==posture & [moveTrajStruct.target]==target).pDPCA.traj(:,1);
        moveTime = moveTrajStruct([moveTrajStruct.posture]==posture & [moveTrajStruct.target]==target).pDPCA.timestamps;        
        holdTraj = holdTrajStruct([holdTrajStruct.posture]==posture & [holdTrajStruct.target]==target).pDPCA.traj(:,1);
        holdTime = holdTrajStruct([holdTrajStruct.posture]==posture & [holdTrajStruct.target]==target).pDPCA.timestamps;
        holdTime = holdTime + timeOffset;        
        if target == plotTargetList(1)
            plot(moveTime,moveTraj,'--','LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
            plot(holdTime,holdTraj,'--','LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
        elseif target == plotTargetList(2)
            plot(moveTime,moveTraj,'LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
            plot(holdTime,holdTraj,'LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
        end     
        ylabel('Posture Dim 1 (a.u.)')
        
        %Hand (Along target axis)
        subplot(2,1,2);hold on
        moveTraj = moveTrajStruct([moveTrajStruct.posture]==posture & [moveTrajStruct.target]==target).targetAxis.traj(:,1);
        moveTime = moveTrajStruct([moveTrajStruct.posture]==posture & [moveTrajStruct.target]==target).targetAxis.timestamps;
        holdTraj = holdTrajStruct([holdTrajStruct.posture]==posture & [holdTrajStruct.target]==target).targetAxis.traj(:,1);
        holdTime = holdTrajStruct([holdTrajStruct.posture]==posture & [holdTrajStruct.target]==target).targetAxis.timestamps;
        holdTime = holdTime + timeOffset;      
        if target == plotTargetList(1)
            plot(moveTime,moveTraj,'--','LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
            plot(holdTime,holdTraj,'--','LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
        elseif target == plotTargetList(2)
            plot(moveTime,moveTraj,'LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
            plot(holdTime,holdTraj,'LineWidth',traceLineWidth,'Color',pcmap(postureInd,:));
        end        
        ylabel('Hand Pos (mm)')
    end
    postureInd = postureInd + 1;
end

subplot(2,1,1)
    xticks([-200, 0, 300, 350, 550])
    xticklabels({})
    ylabel(['Posture',newline,'Dim 1 (a.u.)'])
    ax = gca;
    rectangle('Position',[0,ax.YLim(1),300,ax.YLim(2)-ax.YLim(1)], 'FaceColor',[0.7 0.7 0.7 .25], 'EdgeColor',[1 1 1 0])
    ax.YTick = [ax.YTick(1) ax.YTick(end)];
    ax.TickDir = 'out';
    ax.FontName = 'arial';
    ax.FontSize = traceFs;
    
subplot(2,1,2)
    xticks([-200, 0, 300, 350, 550])
    xticklabels({'Move-200','Move','Move+300','Trial End-200','Trial End'})
    xtickangle(45) 
    xlabel('time (ms)')
    ylabel(['Hand',newline,'Pos (mm)'])
    ax = gca;
    ax.YTick = [ax.YTick(1) ax.YTick(end)];
    ylim = ax.YLim;
    rectangle('Position',[0,ax.YLim(1),300,ax.YLim(2)-ax.YLim(1)], 'FaceColor',[0.7 0.7 0.7 .25], 'EdgeColor',[1 1 1 0])
    ax.YLim = ylim;
    ax.TickDir = 'out';
    ax.TickDir = 'out';
    ax.FontName = 'arial';
    ax.FontSize = traceFs;
    
    if saveFig
        saveas(gcf,fullfile(saveDir,[dataset,'_dPCA_Traces.svg']));
    end
    

    
end