clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 3';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set Parameters
   numPCsToKeep = 10;     %Num PCs to project data into before analysis
   numPdims = 2;
   numTdims = 2;
   inclAvgDot = true; %If true, make traces transparent and add dot for average on top

%% Run loop for each dataset      
    %Load data
    dataset = 'E20200318';
    [Data,zScoreParams] = loadData(dataset);
    [Data] = removeShortBCIandIsoTrials(Data,dataset);
    %Get trajStruct
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);                 
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
    %Project all data down to top PCs
    [allTraj] = collectAllAvgTraj(trajStruct);
    [trajStruct] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep);
    
%% Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(numPostures);      
    
%% Fit P and T spaces, components for each group
    [minNumTimestamps] = getMinNumTimestamps(trajStruct);    
    [pSig,tSig] = getPandTsig(trajStruct,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
    pSigReshape = reshape(squeeze(pSig),[numPostures*minNumTimestamps,numPCsToKeep]);
    tSigReshape = reshape(squeeze(tSig),[numTargets*minNumTimestamps,numPCsToKeep]);
    [pDims,~,~,~,explainedP,pSigMu] = pca(pSigReshape); 
    [tDims,~,~,~,explainedT,tSigMu] = pca(tSigReshape);
 
%% Plot
        fs = 6;
        dotSize = 12;%3;
        figScale = 0.24;
        lineWidth = .65;
        dotLineWidth = .2;
        dotAlpha = .95;
        lineAlpha = .6;
        trajEndDotSize = 2;
        trajEndDotAlpha = .6;
        
        %Plot T signal in T subspace    
        f = figure; 
        f.Position = [200 200 f.Position(3)*figScale f.Position(4)*figScale]; 
        hold on
        targetInd = 1;
        for target = targetList
            traj = tSigReshape((targetInd-1)*minNumTimestamps+1:targetInd*minNumTimestamps,:);
            proj = traj*tDims;
            if inclAvgDot
                trajAvg = mean(proj(:,[1,2]));
                plot(proj(:,1),proj(:,2),'LineWidth',lineWidth,'Color',[tcmap(target,:),lineAlpha]);
                scatter(proj(end,1),proj(end,2),trajEndDotSize,tcmap(target,:),'filled');
                alpha(trajEndDotAlpha);
                scatter(trajAvg(:,1),trajAvg(:,2),dotSize,tcmap(target,:),'filled');
                alpha(dotAlpha);
                %plot(trajAvg(:,1),trajAvg(:,2),'o','MarkerSize',dotSize,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',[0 0 0],'LineWidth',dotLineWidth) 
                %plot(trajAvg(:,1),trajAvg(:,2),'o','MarkerSize',dotSize,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',[0 0 0],'LineWidth',dotLineWidth)    
            else
                plot(proj(:,1),proj(:,2),'LineWidth',lineWidth,'Color',tcmap(target,:));
            end
            targetInd = targetInd + 1;
        end
        %xlabel('Target Dim 1'); ylabel('Target Dim 2')
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        axis equal; 
        set(gca,'TickDir','out');
        ax = gca; TSubXLims = ax.XLim; TSubYLims = ax.YLim;
        ax.XTick = [ax.XTick(1) ax.XTick(end)];
        ax.YTick = [ax.YTick(1) ax.YTick(end)];
        TSubXRange = TSubXLims(2)-TSubXLims(1); TSubYRange = TSubYLims(2)-TSubYLims(1);
        %set(gca,'Layer','bottom')
        
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_TSigTSub.svg']));
        end
                     
        
        %Plot P signal in P subspace
        f = figure; 
        f.Position = [200 200 f.Position(3)*figScale f.Position(4)*figScale]; 
        hold on
        postureInd = 1;
        for posture = postureList
            traj = pSigReshape((postureInd-1)*minNumTimestamps+1:postureInd*minNumTimestamps,:);
            proj = traj*pDims;
            if inclAvgDot
                trajAvg = mean(proj(:,[1,2]));
                plot(-proj(:,1),proj(:,2),'LineWidth',lineWidth,'Color',[pcmap(posture,:),lineAlpha]);
                scatter(-proj(end,1),proj(end,2),trajEndDotSize,pcmap(posture,:),'filled');
                alpha(trajEndDotAlpha);
                scatter(-trajAvg(:,1),trajAvg(:,2),dotSize,pcmap(posture,:),'filled');
                alpha(dotAlpha);
                %plot(trajAvg(:,1),trajAvg(:,2),'o','MarkerSize',dotSize,'MarkerFaceColor',pcmap(posture,:),'MarkerEdgeColor',[0 0 0],'LineWidth',dotLineWidth)
            else
                plot(-proj(:,1),proj(:,2),'LineWidth',lineWidth,'Color',pcmap(posture,:));
            end
            postureInd = postureInd + 1;
        end
        %xlabel('Posture Dim 1'); ylabel('Posture Dim 2')
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        axis equal; 
        set(gca,'TickDir','out');
        ax = gca; PSubXLims = ax.XLim; PSubYLims = ax.YLim;
        ax.XTick = [ax.XTick(1) ax.XTick(end)];
        ax.YTick = [ax.YTick(1) ax.YTick(end)];
        PSubXRange = PSubXLims(2)-PSubXLims(1); PSubYRange = PSubYLims(2)-PSubYLims(1);
        rectangle('Position',[ax.XLim(1),ax.YLim(1),ax.XLim(2)-ax.XLim(1),ax.YLim(2)-ax.YLim(1)], 'FaceColor',[0.7 0.7 0.7 .25], 'EdgeColor',[1 1 1 0])
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PSigPSub.svg']));
        end

        %Plot P signal in T subspace 
        f = figure; 
        f.Position = [200 200 f.Position(3)*figScale f.Position(4)*figScale]; 
        hold on
        postureInd = 1;
        for posture = postureList
            traj = pSigReshape((postureInd-1)*minNumTimestamps+1:postureInd*minNumTimestamps,:);
            proj = traj*tDims;
            if inclAvgDot
                trajAvg = mean(proj(:,[1,2]));
                plot(proj(:,1),proj(:,2),'LineWidth',lineWidth,'Color',[pcmap(posture,:),lineAlpha]);
                scatter(proj(end,1),proj(end,2),trajEndDotSize,pcmap(posture,:),'filled');
                alpha(trajEndDotAlpha);
                scatter(trajAvg(:,1),trajAvg(:,2),dotSize,pcmap(posture,:),'filled');
                alpha(dotAlpha);
                %plot(trajAvg(:,1),trajAvg(:,2),'o','MarkerSize',dotSize,'MarkerFaceColor',pcmap(posture,:),'MarkerEdgeColor',[0 0 0],'LineWidth',dotLineWidth)
            else
                plot(proj(:,1),proj(:,2),'LineWidth',lineWidth,'Color',pcmap(posture,:));
            end
            postureInd = postureInd + 1;
        end
        %xlabel('Target Dim 1'); ylabel('Target Dim 2')
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        ax = gca; xlim = ax.XLim; ylim = ax.YLim;
        xmid = mean(xlim); ymid = mean(ylim);
        ax.XLim = [xmid-PSubXRange/2 xmid+PSubXRange/2]; 
        ax.YLim = [ymid-PSubYRange/2 ymid+PSubYRange/2]; 
        ax.XTick = [ax.XTick(1) ax.XTick(end)];
        ax.YTick = [ax.YTick(1) ax.YTick(end)];
        set(gca,'TickDir','out');
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PSigTSub.svg']));
        end

        
        %Plot T signal in P subspace 
        f = figure; 
        f.Position = [200 200 f.Position(3)*figScale f.Position(4)*figScale]; 
        hold on
        targetInd = 1;
        for target = targetList
            traj =  tSigReshape((targetInd-1)*minNumTimestamps+1:targetInd*minNumTimestamps,:);
            proj = traj*pDims;        
            if inclAvgDot
                trajAvg = mean(proj(:,[1,2]));
                plot(-proj(:,1),proj(:,2),'LineWidth',lineWidth,'Color',[tcmap(target,:),lineAlpha]);
                scatter(-proj(end,1),proj(end,2),trajEndDotSize,tcmap(target,:),'filled');
                alpha(trajEndDotAlpha);
                scatter(-trajAvg(:,1),trajAvg(:,2),dotSize,tcmap(target,:),'filled');
                alpha(dotAlpha);
                %plot(trajAvg(:,1),trajAvg(:,2),'o','MarkerSize',dotSize,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',[0 0 0],'LineWidth',dotLineWidth)
            else
                plot(-proj(:,1),proj(:,2),'LineWidth',lineWidth,'Color',tcmap(target,:));
            end          
            targetInd = targetInd + 1;
        end
        %xlabel('Posture Dim 1'); ylabel('Posture Dim 2')
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        ax = gca; xlim = ax.XLim; ylim = ax.YLim;

        xmid = mean(xlim); ymid = mean(ylim);
        ax.XLim = [xmid-TSubXRange/2 xmid+TSubXRange/2]; 
        ax.YLim = [ymid-TSubYRange/2 ymid+TSubYRange/2]; 
        ax.XTick = [ax.XTick(1) ax.XTick(end)];
        ax.YTick = [ax.YTick(1) ax.YTick(end)];
        set(gca,'TickDir','out');
        rectangle('Position',[ax.XLim(1),ax.YLim(1),ax.XLim(2)-ax.XLim(1),ax.YLim(2)-ax.YLim(1)], 'FaceColor',[0.7 0.7 0.7 .25], 'EdgeColor',[1 1 1 0])
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_TSigPSub.svg']));
        end
