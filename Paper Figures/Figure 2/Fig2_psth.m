clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 2';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Load data
    dataset = 'E20200316';
    [Data,zScoreParams] = loadData(dataset);
    [Data] = removeShortBCIandIsoTrials(Data,dataset);
    
%% Get trajStruct
    %Get trajStruct - adjust for (1) time window (-50 to +250ms rel. Go) & (2) no z-score 
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    trajFields = {'smoothFR'};
    trialInclStates(1).inclStates = {{'state','Step 1','first',-50},{'state','Step 1','first',250}};    
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',true);    
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 

%% Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(numPostures);  

%% Load anova table for dataset
    anovaResultStruct = load(fullfile('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S1 - single units',dataset,'anovaResult.mat'));
    anovaResultStruct = anovaResultStruct.resultStruct;
    
%% Plot all PSTHs
    chGroups = {1:49,50:87};       
    %Plot all targets for posture 1
    posture = 1;
    for chGroup = 1:2
        chList = chGroups{chGroup};
        f = figure; fs = 14; f.Position = [200,200,1225,464];        
        for target = targetList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
            chInd = 1;
            for ch = chList
                subplot(7,7,chInd);
                plot(time,traj(:,ch),'LineWidth',1.5,'Color',tcmap(target,:));
                hold on;
                if target == targetList(end)
                    ax = gca;
                    curYLims = ax.YLim;        
                    text(0,curYLims(2),num2str(ch))
                end
                chInd = chInd + 1;
            end
        end
    end
     
    %Plot all postures for target 1
    target = 1;
    for chGroup = 1:2
        chList = chGroups{chGroup};
        f = figure; fs = 14; f.Position = [200,200,1225,464];  
        for posture = postureList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
            chInd = 1;
            for ch = chList
                subplot(7,7,chInd);
                plot(time,traj(:,ch),'LineWidth',1.5,'Color',pcmap(posture,:));
                hold on;
                if posture == postureList(end)
                    ax = gca;
                    curYLims = ax.YLim;        
                    text(0,curYLims(2),num2str(ch))
                end
                chInd = chInd + 1;
            end
        end
    end
    
%% Select units for Fig 2
    %Plot a select few neurons
    %postureAndTarget = [3,42];  
    %targetOnly = [29,34,41,43,46];
    %postureOnly = [49,72,70];
    %prettiest = [3,34,41,42,43,82];
    mixedUnits = find([anovaResultStruct.tuning]==3);
    targetUnits = find([anovaResultStruct.tuning]==1);
    postureUnits = find([anovaResultStruct.tuning]==2);
    chList = [3,34,49];
    
%% Plot selected units 
%chList = [5,46,87];
    %Plot all postures for target 1
    f = figure; f.Position = [100,-200,518,950];
    target = 1;
    for posture = postureList
        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
        CI = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.CI95;
        SEM = CI./1.96;
        time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
        chInd = 1;
        for ch = chList
            subplot(3,1,chInd)
            shadedErrorBar(time,traj(:,ch),SEM(:,ch),'lineprops',{'LineWidth',5,'Color',pcmap(posture,:)});
            hold on;
            chInd = chInd + 1;
        end
    end
    format3NeuronPSTH();
    if saveFig
        saveas(gcf,fullfile(saveDir,[dataset,'_posture_PSTH.svg']));
    end
    
    %Plot all targets for posture 1
    f = figure; f.Position = [100,-200,518,950];
    posture = 1;
    for target = targetList
        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
        CI = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.CI95;
        SEM = CI./1.96;
        time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
        chInd = 1;
        for ch = chList
            subplot(3,1,chInd)
            shadedErrorBar(time,traj(:,ch),SEM(:,ch),'lineprops',{'LineWidth',5,'Color',tcmap(target,:)});
            hold on;
            chInd = chInd + 1;
        end
    end
    format3NeuronPSTH();
    if saveFig
        saveas(gcf,fullfile(saveDir,[dataset,'_target_PSTH.svg']));
    end
    
%% Local functions
    function [] = format3NeuronPSTH()
        for chInd = 1:3
            subplot(3,1,chInd); hold on
            ax = gca; xlimits = ax.XLim; ylimits = ax.YLim;
            ax.TickDir = 'out';
            set(gca,'fontname','arial'); set(gca,'fontsize',14);
            xticklabels({});
            if chInd == 3
                xticks([0,100,200])
                xlabel('time (ms)')
            end 
        end
    end