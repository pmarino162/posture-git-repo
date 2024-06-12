clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S7 - posture dim during reaching';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Parameters
    numPCsToKeep = 10;
    regularize = true;
    
%% Run loop for each dataset  
reachDatasetList = {'N20190226','R20200221','E20210706'};

for datasetList = {'E20210706'}
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
    % Keep postures 1-4 only for earl reaching
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
    pDPCA = W(:,postureDims);
    
    
    %% Get marker axis 
    switch dataset
        case {'N20190226','R20200221'}
            targetID = [Data.targetData];
            targetID = [targetID.targetID];

            target4Loc = Data(find(targetID==4,1,'first')).targetData.targetLoc(1,1:2);
            target8Loc = Data(find(targetID==8,1,'first')).targetData.targetLoc(1,1:2);
            targetAxis = (target4Loc-target8Loc)./norm(target4Loc-target8Loc);
        case{'E20210706'}
            targetAxis = [1 0];
            targetAxis = [-0.7071 0.7071];
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
        moveTrajStruct(i).targetAxis.traj = moveTrajStruct(i).avgMarker.traj(:,1:2)*targetAxis';
        moveTrajStruct(i).targetAxis.timestamps = moveTrajStruct(i).avgMarker.timestamps;
             
        holdTrajStruct(i).avgPCA.traj = (holdTrajStruct(i).avgZSmoothFR.traj-Mu)*PCs;
        holdTrajStruct(i).pPCA.traj = holdTrajStruct(i).avgPCA.traj(:,1:numPCsToKeep)*pDims;
        holdTrajStruct(i).pPCA.timestamps = holdTrajStruct(i).avgZSmoothFR.timestamps;     
        holdTrajStruct(i).pDPCA.traj = holdTrajStruct(i).avgZSmoothFR.traj*pDPCA;
        holdTrajStruct(i).pDPCA.timestamps = holdTrajStruct(i).avgZSmoothFR.timestamps;      
        holdTrajStruct(i).targetAxis.traj = holdTrajStruct(i).avgMarker.traj(:,1:2)*targetAxis';
        holdTrajStruct(i).targetAxis.timestamps = holdTrajStruct(i).avgMarker.timestamps;
        
        %Add all traces to struct
        for j = 1:size(trajStruct(i).allZSmoothFR,2)
            moveTrajStruct(i).allPCA(j).traj = (moveTrajStruct(i).allZSmoothFR(j).traj-Mu)*PCs;
            moveTrajStruct(i).allPPCA(j).traj = moveTrajStruct(i).allPCA(j).traj(:,1:numPCsToKeep)*pDims;
            moveTrajStruct(i).allPPCA(j).timestamps = moveTrajStruct(i).allZSmoothFR(j).timestamps;
            moveTrajStruct(i).allTargetAxis(j).traj = moveTrajStruct(i).allMarker(j).traj(:,1:2)*targetAxis';
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
    
%% Plot correlation
    plotInd = 1;
    for plotPostureList = [2,7]
        for plotTargetList = [4,8]
    %         plotPostureList = 1;
    %         plotTargetList = 8;
            allData = zeros(1,2);
            allDataInd = 1;
            f = figure; hold on
            f.Position = [200 200 f.Position(3)*corrFigScale f.Position(4)*corrFigScale];
            for posture = plotPostureList
                for target = plotTargetList
                    trajStructInd = find([moveTrajStruct.posture]==posture & [moveTrajStruct.target]==target);
                    for i = 1:size(trajStruct(trajStructInd).allZSmoothFR,2)
                        time = moveTrajStruct(trajStructInd).allPPCA(i).timestamps;
                        timeMask = time >= 0;
                        neuralTraj = moveTrajStruct(trajStructInd).allPPCA(i).traj(timeMask,1);
                        handTraj = moveTrajStruct(trajStructInd).allTargetAxis(i).traj(timeMask,1);
                        plot(neuralTraj,handTraj,'.k','MarkerSize',corrDotSize);
                        numPts = sum(timeMask);
                        allData(allDataInd:allDataInd+numPts-1,1) = neuralTraj;
                        allData(allDataInd:allDataInd+numPts-1,2) = handTraj;
                        allDataInd = allDataInd + numPts;
                    end
                end
            end
            xlabel('Posture Dim 1 (a.u.)')
            ylabel('Hand Pos (mm)')

            [R,P] = corrcoef(allData,'Alpha',0.05);
            ax = gca;
            xlim = ax.XLim;
            ylim = ax.YLim;
            xrange = xlim(2)-xlim(1);
            xmid = xrange/2;
            yrange = ylim(2)-ylim(1);
            ymid = yrange/2;
             ax.XTick = [ax.XTick(1) ax.XTick(end)];
              ax.YTick = [ax.YTick(1) ax.YTick(end)];
            text(xlim(1)+xrange*.8,ylim(1)+yrange*.8,['\rho = ',num2str(R(1,2)),newline,'p = ',num2str(P(1,2))],'FontSize',corrFs)
            ax.FontName = 'arial';
            ax.FontSize = corrFs;


            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'ind_',num2str(plotInd),'_PCA_Corr.svg']));
            end

            plotInd = plotInd + 1;
        end
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
    
%% Plot correlation
    plotInd = 1;
    for plotPostureList = [2,7]
        for plotTargetList = [4,8]
    allData = zeros(1,2);
    allDataInd = 1;
    f = figure; hold on
    f.Position = [200 200 f.Position(3)*corrFigScale f.Position(4)*corrFigScale];
    for posture = plotPostureList
        for target = plotTargetList
            trajStructInd = find([moveTrajStruct.posture]==posture & [moveTrajStruct.target]==target);
            for i = 1:size(trajStruct(trajStructInd).allZSmoothFR,2)
                time = moveTrajStruct(trajStructInd).allPPCA(i).timestamps;
                timeMask = time >= 0;
                neuralTraj = moveTrajStruct(trajStructInd).allDPCA(i).traj(timeMask,1);
                handTraj = moveTrajStruct(trajStructInd).allTargetAxis(i).traj(timeMask,1);
                plot(neuralTraj,handTraj,'.k','MarkerSize',corrDotSize);
                numPts = sum(timeMask);
                allData(allDataInd:allDataInd+numPts-1,1) = neuralTraj;
                allData(allDataInd:allDataInd+numPts-1,2) = handTraj;
                allDataInd = allDataInd + numPts;
            end
        end
    end
    xlabel('Posture Dim 1 (a.u.)')
    ylabel('Hand Pos (mm)')

    [R,P] = corrcoef(allData,'Alpha',0.05);
    ax = gca;
    xlim = ax.XLim;
    ylim = ax.YLim;
    xrange = xlim(2)-xlim(1);
    xmid = xrange/2;
    yrange = ylim(2)-ylim(1);
    ymid = yrange/2;
     ax.XTick = [ax.XTick(1) ax.XTick(end)];
      ax.YTick = [ax.YTick(1) ax.YTick(end)];
    text(xlim(1)+xrange*.8,ylim(1)+yrange*.8,['\rho = ',num2str(R(1,2)),newline,'p = ',num2str(P(1,2))],'FontSize',corrFs)
    ax.FontName = 'arial';
    ax.FontSize = corrFs;
    
    if saveFig
        saveas(gcf,fullfile(saveDir,[dataset,'ind_',num2str(plotInd),'_dPCA_Corr.svg']));
    end
    plotInd = plotInd + 1;
        end
    end
    
    
end