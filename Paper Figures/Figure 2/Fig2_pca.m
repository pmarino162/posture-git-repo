clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 2';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Parameters
    numPCsToKeep = 10;
    
%% Load data
    dataset = 'E20200318';
    [Data,zScoreParams] = loadData(dataset);
    [Data] = removeShortBCIandIsoTrials(Data,dataset);
    
%% Get trajStruct - fit PCA on usual anlaysis window, project extra data into these dims for fig S2
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    %trajStruct - for Fig 2 and PC fitting
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
    %exTrajStruct - for Fig S2 projection
    trialInclStates(1).inclStates = {{'state','Step 1','first',-50},{'state','Step 1','first',250}}; 
    exTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);

%% Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(numPostures);   
    
%% Fit PCs
    [allTraj,Mu] = collectAllAvgTraj(trajStruct);
    [trajStruct,PCs] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep);
    [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
    
%% Add Projections to exTrajStruct
    [allExTraj,ExMu] = collectAllAvgTraj(exTrajStruct);
    totalVar = trace(cov(allExTraj));
    for i = 1:size(exTrajStruct,2)
        %Add average traces to exTrajStruct
        exTrajStruct(i).avgPCA.traj = (exTrajStruct(i).avgZSmoothFR.traj-Mu)*PCs(:,1:numPCsToKeep);
        exTrajStruct(i).avgPCA.VAF =  100.*(diag(cov(allExTraj*PCs(:,1:numPCsToKeep)))')./totalVar;
        exTrajStruct(i).avgPCA.timestamps = exTrajStruct(i).avgZSmoothFR.timestamps;
        %Compute CI, add to exTrajStruct
        %adjust CI so that it uses available number of trials at each time point
        numTrials = size(exTrajStruct(i).allZSmoothFR,2);
        maxNumTimestamps = 0;
        for trial = 1:numTrials
           numTimestamps = size(exTrajStruct(i).allZSmoothFR(trial).traj,1);
           if numTimestamps > maxNumTimestamps
               maxNumTimestamps = numTimestamps;
           end
        end
        tempAllTraj = NaN(maxNumTimestamps,numTrials,numPCsToKeep);
        for trial = 1:numTrials
            traj = exTrajStruct(i).allZSmoothFR(trial).traj;
            numTimestamps = size(traj,1);
            tempAllTraj(1:numTimestamps,trial,:) = (traj-Mu)*PCs(:,1:numPCsToKeep);
        end
        for timestamp = 1:maxNumTimestamps
            timestampMat = squeeze(tempAllTraj(timestamp,:,:)); %All trials for all channels for current timestamp
            tempNumTrials = sum(~isnan(timestampMat(:,1)));
            exTrajStruct(i).avgPCA.CI(timestamp,:) = 1.96.*nanstd(timestampMat)./sqrt(tempNumTrials);
        end
    end
        
%% Plot
    %Orthographic for Fig 2
    plotPostureList = 1:5;
    plotTargetList = 1:8;
    fs = 14;
    figure
    %f = figure; f.Position = [50 50 200 200];
    hold on
    timePts = 1:minNumTimestamps;
    xDim = 1; yDim = 2; zDim = 3;
    for posture = plotPostureList
        for target = plotTargetList
            if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPCA.traj(timePts,:); 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end
    grid on
    xticklabels({}); yticklabels({}); zticklabels({}); 
    VAF = round(trajStruct(1).avgPCA.VAF);
    xlabel(['PC',num2str(xDim),' (',num2str(VAF(xDim)),'%)'])
    ylabel(['PC',num2str(yDim),' (',num2str(VAF(yDim)),'%)'])
    zlabel(['PC',num2str(zDim),' (',num2str(VAF(zDim)),'%)'])
    view([20 10])
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
            
    %Save
    if saveFig
        saveas(gcf,fullfile(saveDir,[dataset,'_PC_orth.fig']));
    end
    
    
%% Traces for Fig S2
    plotPostureList = [1,3,5];
    plotTargetList = [1:8];
    fs = 4;
    for target = plotTargetList
        f = figure; f.Position = [50 50 200 80];
        for posture = plotPostureList
            traj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgPCA.traj;
            CI = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgPCA.CI;
            time = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgPCA.timestamps;  
            for dim = 1:6
                subplot(2,3,dim); hold on;
                shadedErrorBar(time,traj(:,dim),CI(1:length(time),dim),'lineprops',{'LineWidth',1,'Color',pcmap(posture,:)});
            end
        end    
        
        %Format
        maxYRange = 0;
        for dim = 1:6
            subplot(2,3,dim)
            ax = gca;
            ylimits = ax.YLim;
            yRange = ylimits(2)-ylimits(1);
            if yRange > maxYRange
                maxYRange = yRange;
            end
        end
        for dim = 1:6
            subplot(2,3,dim);
            ax = gca;
            ylimits = ax.YLim;
            yMid = ylimits(1) + 0.5*(ylimits(2)-ylimits(1));
            ylim([yMid-maxYRange/2 yMid+maxYRange/2])
            
            plot([0 0],[yMid-maxYRange/2 yMid+maxYRange/2],':k','LineWidth',1)
            ylabel(['PC ',num2str(dim)])
            yticks([min(ax.YTick) max(ax.YTick)]);
            ax.TickDir = 'out';
            xticks([-50 0 250]);
            if dim < 4
                xticklabels({})
            else
                xticklabels({'-50','','250'});
            end
            if dim == 5
                %xlabel('time (ms)')
            end
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
        end           
     
        %Save
        if saveFig
            saveas(gcf,fullfile(saveDir,['Figure S2\',dataset,'_target_',num2str(target),'_PC_traces.svg']));
        end
    end
