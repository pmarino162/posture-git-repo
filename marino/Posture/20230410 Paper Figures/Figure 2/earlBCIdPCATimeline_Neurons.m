clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20230203\orthographic dpca pop traj';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
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
    exTrajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    

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
     plotTargetList = [3,7];
    %All Neurons, numbered
    for channelGroup = 1
        if channelGroup == 1
            dimList = 1:49;
        else
            dimList = 50:87;
        end

        f = figure; fs = 14;
        f.Position = [200,200,1225,464];
        postureInd = 1;
        for posture = postureList
            for target = plotTargetList
                exTraj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.traj(1:numPts,:);
                exCI = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.CI95(1:numPts,:);
                exCI = exCI./1.96;
                exTime = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.timestamps(1:numPts);

                dimInd = 1;
                for dim = dimList
                    subplot(7,7,dimInd); hold on
                    if target == 1 || target == 3
                        plot(exTime,exTraj(:,dim),'LineWidth',1.5,'Color',pcmap(posture,:));
                        hold on;
                    elseif target == 5 || target == 7
                        plot(exTime,exTraj(:,dim),'--','LineWidth',1.5,'Color',pcmap(posture,:));
                        hold on;
                    end
                    if posture == postureList(end) && target == plotTargetList(end)
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
    f = figure; fs = 14;
    f.Position = [100,-200,518,800];
    postureInd = 1;
    dimList = [3,49,78,76];
    for posture = postureList
        for target = plotTargetList

            exTraj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.traj(1:numPts,:);
            exCI = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.CI95(1:numPts,:);
            exCI = exCI./1.96;
            exTime = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.timestamps(1:numPts);

            dimInd = 1;
            for dim = dimList
                subplot(4,1,dimInd); hold on
                if target == 1 || target == 3
                    shadedErrorBar(exTime,exTraj(:,dim),exCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                    hold on;
                elseif target == 5 || target == 7
                    shadedErrorBar(exTime,exTraj(:,dim),exCI(:,dim),'lineprops',{'--','LineWidth',1.5,'Color',pcmap(posture,:)});
                    hold on; 
                end
                dimInd = dimInd + 1;
            end
        end
        postureInd = postureInd + 1;
    end

    dimInd = 1;
    for dim = dimList(1:4)
        subplot(4,1,dimInd); hold on
        ax = gca; xlimits = ax.XLim; ylimits = ax.YLim;
        ax.TickDir = 'out';
        set(gca,'fontname','arial'); set(gca,'fontsize',fs);
        if ismember(dimInd,[1,2,3])
            xticklabels({});
        end
        if dimInd == 4
            xticks([0,100,200])
            xlabel('time (ms)')
        end 

        dimInd = dimInd + 1;
    end

  