clear; clc; clf; close all;
    
%% Setup save dir   
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20210930 - first model\';

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = rainbow2d;
    
%% Load Data 
    binWidth = 25;
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
    exclCh =  [44 87 88 77 78 71 67 69 118];
    getSorts = false;
    Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'binWidth',binWidth,'exclZero',false,'getMarker',true,'centerMarker',false,'getKin',true);
    [Data] = cleanData20210706(Data);
    Data = Data([Data.trialStatus]==1);
    [Data,postureIDs] = labelPostures20210706(Data);
    allTrialPostures = [Data.conditionData];
    allTrialPostures = [allTrialPostures.postureID];
%     keepPostures = [1,2,4,5];
    keepPostures = [1:7];
    Data = Data(ismember(allTrialPostures,keepPostures));
    
%% Plot target posture map
    conditionData = [Data.conditionData];
    targetData = [Data.targetData];
    figure
    hold on
    for posture = [1,4,7]
        for target = 1:8
            tempData = Data([conditionData.postureID]==posture & [targetData.targetID]==target);
            if ~isempty(tempData)
                workspaceCenter = tempData(1).targetData.workspaceCenter;
                targetLoc = tempData(1).targetData.targetLoc;
                targetLoc = targetLoc + workspaceCenter;
                plot(targetLoc(1),targetLoc(2),'.','MarkerSize',20,'Color',tcmap(target,:));
                if target == 1
                   text(workspaceCenter(1),workspaceCenter(2),num2str(posture)) 
                end
            end
        end
    end
    axis equal
    xlabel('x (mm)')
    ylabel('y (mm)')

%% Plot reaction times by condition

    figure
    hold on
    for posture = keepPostures
        for target = 1:8
            trialNum = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).trialNum];
            kinData = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).kinData];
            if ~isempty(kinData)
                rxnTimes = [kinData.reactTimeMoveOnset];
                jitter = 0.25*rand(1,length(rxnTimes));
                plot(rxnTimes,posture*ones(1,length(rxnTimes))+jitter,'.','MarkerSize',7,'Color',tcmap(target,:));
%                 text(rxnTimes,posture*ones(1,length(rxnTimes))+jitter,num2str(trialNum'));
            end
        end
    end
    xlabel('Reaction Time (ms)')
    ylabel('Posture Label')
    
%% Condition median reaction times
    figure
    hold on
    for posture = keepPostures
        for target = 1:8
            trialNum = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).trialNum];
            kinData = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).kinData];
            if ~isempty(kinData)
                rxnTimes = median([kinData.reactTimeMoveOnset]);
                jitter = 0.1*rand(1,1);
                plot(rxnTimes,posture+jitter,'.','MarkerSize',20,'Color',tcmap(target,:));
%                 text(reachTimes,posture*ones(1,length(reachTimes))+jitter,num2str(trialNum'))
            end
        end
    end
    xlabel('Median Reaction Time (ms)')
    ylabel('Posture Label')
    
%% Plot reaction time histograms by posture
    %Get histograms
    histStruct = struct('posture',[],'target',[],'hX',[],'hY',[],'totalTrials',[]);
    structInd = 1;
    for posture = keepPostures
        for target = 1:8
            kinData = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).kinData];
            if ~isempty(kinData)
                rxnTimes = [kinData.reactTimeMoveOnset];
                h = histfit(rxnTimes);
                set(gcf,'Visible','off');
                hX = h(2).XData;
                hY = h(2).YData;
                histStruct(structInd).posture = posture;
                histStruct(structInd).target = target;
                histStruct(structInd).hX = hX;
                histStruct(structInd).hY = hY;
                histStruct(structInd).totalTrials = size(kinData,2);
                structInd = structInd + 1;
            end
        end
    end
    
    figure
    hold on
    for posture = keepPostures
        for target = 1:8
            tempData = histStruct([histStruct.posture]==posture & [histStruct.target]==target);
            if ~isempty(tempData)
                totalTrials = tempData.totalTrials;
                hX = tempData.hX;
                hY = tempData.hY./totalTrials;
                plot(hX,posture+hY,'LineWidth',3,'Color',tcmap(target,:));
%                 text(reachTimes,posture*ones(1,length(reachTimes))+jitter,num2str(trialNum'))
            end
        end
    end
    xlabel('Reaction Time (ms)')
    ylabel('Posture Label')



%% Plot reach times by condition
    figure
    hold on
    for posture = keepPostures
        for target = 1:8
            trialNum = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).trialNum];
            kinData = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).kinData];
            if ~isempty(kinData)
                reachTimes = [kinData.reachTime];
                jitter = 0.25*rand(1,length(reachTimes));
                plot(reachTimes,posture*ones(1,length(reachTimes))+jitter,'.','MarkerSize',7,'Color',tcmap(target,:));
%                 text(reachTimes,posture*ones(1,length(reachTimes))+jitter,num2str(trialNum'))
            end
        end
    end
    xlabel('Reach Time (ms)')
    ylabel('Posture Label')
    
    
%% Condition median reach times
    figure
    hold on
    for posture = keepPostures
        for target = 1:8
            trialNum = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).trialNum];
            kinData = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).kinData];
            if ~isempty(kinData)
                reachTime = median([kinData.reachTime]);
                jitter = 0.1*rand(1,1);
                plot(reachTime,posture+jitter,'.','MarkerSize',20,'Color',tcmap(target,:));
%                 text(reachTimes,posture*ones(1,length(reachTimes))+jitter,num2str(trialNum'))
            end
        end
    end
    xlabel('Median Reach Time (ms)')
    ylabel('Posture Label')
    
    
    %% Plot reach time histograms by posture
    %Get histograms
    histStruct = struct('posture',[],'target',[],'hX',[],'hY',[],'totalTrials',[]);
    structInd = 1;
    for posture = keepPostures
        for target = 1:8
            kinData = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).kinData];
            if ~isempty(kinData)
                reachTimes = [kinData.reachTime];
                h = histfit(reachTimes);
                set(gcf,'Visible','off');
                hX = h(2).XData;
                hY = h(2).YData;
                histStruct(structInd).posture = posture;
                histStruct(structInd).target = target;
                histStruct(structInd).hX = hX;
                histStruct(structInd).hY = hY;
                histStruct(structInd).totalTrials = size(kinData,2);
                structInd = structInd + 1;
            end
        end
    end
    
    figure
    hold on
    for posture = keepPostures
        for target = 1:8
            tempData = histStruct([histStruct.posture]==posture & [histStruct.target]==target);
            if ~isempty(tempData)
                totalTrials = tempData.totalTrials;
                hX = tempData.hX;
                hY = tempData.hY./totalTrials;
                plot(hX,posture+hY,'LineWidth',2,'Color',tcmap(target,:));
%                 text(reachTimes,posture*ones(1,length(reachTimes))+jitter,num2str(trialNum'))
            end
        end
    end
    xlabel('Reach Time (ms)')
    ylabel('Posture Label')

    
%% Condition median peak speeds (color by posture)
    figure
    hold on
    for target = 1:8
        for posture = keepPostures
            trialNum = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).trialNum];
            kinData = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).kinData];
            if ~isempty(kinData)
                peakSpeeds = median([kinData.peakSpeed]);
                jitter = 0.1*rand(1,1);
                plot(peakSpeeds,target+jitter,'.','MarkerSize',20,'Color',tcmap(posture,:));
%                 text(reachTimes,posture*ones(1,length(reachTimes))+jitter,num2str(trialNum'))
            end
        end
    end
    xlabel('Median Peak Speed (mm/ms)')
    ylabel('Target')
    
    
%% Condition median reach times (color by posture)
    figure
    hold on
    for target = 1:8
        for posture = keepPostures
            trialNum = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).trialNum];
            kinData = [Data([conditionData.postureID]==posture & [targetData.targetID]==target).kinData];
            if ~isempty(kinData)
                reachTime = median([kinData.reachTime]);
                jitter = 0.1*rand(1,1);
                plot(reachTime,target+jitter,'.','MarkerSize',20,'Color',tcmap(posture,:));
%                 text(reachTimes,posture*ones(1,length(reachTimes))+jitter,num2str(trialNum'))
            end
        end
    end
    xlabel('Median Reach Time Time (ms)')
    ylabel('Target')