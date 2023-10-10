clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20230203\orthographic dpca pop traj';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    [pcmap,tcmap,rainbow] = getColorMaps();
    
%% Get trajStruct
    %Load data
    dataset = 'E20200316';
    [Data,zScoreParams] = loadData(dataset);
    
    %Get trajStruct
    binWidth = 25;
    kernelStdDev = 25;
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
    
    %Execution
    trialInclStates(1).inclStates = {{'state','Step 1','first',-50},{'state','Step 2','first',250}};
    exTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    

%% Get numPts
    %Manually enter number of points to consider (execution)
   	numPts = 13;
           
    %Get posture and target lists
    postureList = unique([exTrajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([exTrajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(exTrajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(exTrajStruct,2);

%% Plot
    %Plot all postures for every neuron, then all targets
    plotTargetList = [1:8];
    for channelGroup = 1:2
        if channelGroup == 1
            dimList = 1:49;
        else
            dimList = 50:87;
        end

        f = figure; fs = 14;
        f.Position = [200,200,1225,464];
        postureInd = 1;
        for posture = 1
            for target = [1:8]
                exTraj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.traj(1:numPts,:);
                exCI = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.CI95(1:numPts,:);
                exCI = exCI./1.96;
                exTime = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.timestamps(1:numPts);

                dimInd = 1;
                for dim = dimList
                    subplot(7,7,dimInd); hold on
                    plot(exTime,exTraj(:,dim),'LineWidth',1.5,'Color',tcmap(target,:));
                    hold on;
                    if target == plotTargetList(end)
                        ax = gca;
                        curYLims = ax.YLim;        
                        text(0,curYLims(2),num2str(dim))
                    end
                    dimInd = dimInd + 1;
                end
            end
            postureInd = postureInd + 1;
        end
    end
     
    %All Neurons, numbered
    for channelGroup = 1:2
        if channelGroup == 1
            dimList = 1:49;
        else
            dimList = 50:87;
        end

        f = figure; fs = 14;
        f.Position = [200,200,1225,464];
        postureInd = 1;
        for posture = postureList
            for target = 1
                exTraj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.traj(1:numPts,:);
                exCI = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.CI95(1:numPts,:);
                exCI = exCI./1.96;
                exTime = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.timestamps(1:numPts);

                dimInd = 1;
                for dim = dimList
                    subplot(7,7,dimInd); hold on
                        plot(exTime,exTraj(:,dim),'LineWidth',1.5,'Color',pcmap(posture,:));
                        hold on;
                    if posture == postureList(end) 
                        ax = gca;
                        curYLims = ax.YLim;        
                        text(0,curYLims(2),num2str(dim))
                    end
                    dimInd = dimInd + 1;
                end
            end
            postureInd = postureInd + 1;
        end
    end
    
    
    
    
    
    

    
    %Plot a select few neurons
    postureAndTarget = [3,42];
        targetOnly = [29,34,41,43,46];
        postureOnly = [49,72,70];
      dimList = [82,34,70];
    prettiest = [3,34,41,42,43,82];
    
    dimList = [3,34,49];
    
    dimList = [3,34,65];
    
    f = figure; fs = 14;
     f.Position = [100,-200,518,950];
    postureInd = 1;
    
    plotPostureList = 1:5;
    plotTargetList = 1;
    
    for posture = postureList
        for target = plotTargetList

            exTraj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.traj(1:numPts,:);
            exCI = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.CI95(1:numPts,:);
            exCI = exCI./1.96;
            exTime = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.timestamps(1:numPts);

            dimInd = 1;
            for dim = dimList
                subplot(3,1,dimInd); hold on
                    shadedErrorBar(exTime,exTraj(:,dim),exCI(:,dim),'lineprops',{'LineWidth',5,'Color',pcmap(posture,:)});
                    hold on;
                dimInd = dimInd + 1;
            end
        end
        postureInd = postureInd + 1;
    end

    dimInd = 1;
    for dim = dimList(1:3)
        subplot(3,1,dimInd); hold on
        ax = gca; xlimits = ax.XLim; ylimits = ax.YLim;
        ax.TickDir = 'out';
        set(gca,'fontname','arial'); set(gca,'fontsize',fs);
        if ismember(dimInd,[1,2,3])
            xticklabels({});
        end
        if dimInd == 3
            xticks([0,100,200])
            xlabel('time (ms)')
        end 

        dimInd = dimInd + 1;
    end

   %Plot a select few neurons
   plotPostureList = 1;
   plotTargetList = 1:8;
   
    f = figure; fs = 14;
    f.Position = [100,-200,518,950];
    postureInd = 1;
    for posture = plotPostureList
        for target = plotTargetList

            exTraj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.traj(1:numPts,:);
            exCI = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.CI95(1:numPts,:);
            exCI = exCI./1.96;
            exTime = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.timestamps(1:numPts);

            dimInd = 1;
            for dim = dimList
                subplot(3,1,dimInd); hold on
                shadedErrorBar(exTime,exTraj(:,dim),exCI(:,dim),'lineprops',{'LineWidth',5,'Color',tcmap(target,:)});
                hold on;

                dimInd = dimInd + 1;
            end
        end
        postureInd = postureInd + 1;
    end

    dimInd = 1;
    for dim = dimList(1:3)
        subplot(3,1,dimInd); hold on
        ax = gca; xlimits = ax.XLim; ylimits = ax.YLim;
        ax.TickDir = 'out';
        set(gca,'fontname','arial'); set(gca,'fontsize',fs);
        if ismember(dimInd,[1,2,3])
            xticklabels({});
        end
        if dimInd == 3
            xticks([0,100,200])
            xlabel('time (ms)')
        end 

        dimInd = dimInd + 1;
    end