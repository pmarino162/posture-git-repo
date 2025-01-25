clear; clc; clf; close all;

%% Setup Save
%     savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Task Analysis\Grid Reaching';
    savePath = 'D:\Animals\Earl\2021\06\20210630\Daily Plots\';
    saveFig = false;
    
%% Load Data
    dateStr = '20210706';
%     load('D:\Animals\Earl\2021\06\20210630\Earl20210630_gridReaching_SI_translated.mat')
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat')

%     load('D:\Animals\Earl\2021\06\20210626\Earl20210626_gridReaching_SI_translated.mat');
%     load('D:\Animals\Earl\2021\06\20210628\Earl20210628_gridReaching_SI_translated.mat')
    Data = getDataStruct(Data,'getMarker',true,'getKin',true);

%% Label Postures
    %Get list of targets for grid reaching trials
    gridReachData = Data(cellfun(@(x) strcmpi(x,'GridReaching'),{Data.trialName}));
    targetData = [gridReachData.targetData];
    workspaceCenter = cell2mat({targetData.workspaceCenter}');
    uniWorkspaceCenter = unique(workspaceCenter,'rows');
    sortedYs = sort(unique(uniWorkspaceCenter(:,2)),'descend')';
    sortedXs = sort(unique(uniWorkspaceCenter(:,1)))';
    numCol = size(sortedXs,2);
    numTrials = size(gridReachData,2);
    %Label each trial
    for trial = 1:numTrials
       workspaceCenter = gridReachData(trial).targetData.workspaceCenter(1:2);
       row = find(sortedYs==workspaceCenter(2));
       col = find(sortedXs==workspaceCenter(1));
       gridReachData(trial).conditionData.postureID = (row-1)*numCol + col;
    end
    %Create postureIDs list
    postureIDs = struct('postureID',[],'workspaceCenter',[]);
    for i = 1:size(uniWorkspaceCenter,1)
       workspaceCenter = uniWorkspaceCenter(i,:);
       row = find(sortedYs==workspaceCenter(2));
       col = find(sortedXs==workspaceCenter(1));
       postureID = (row-1)*numCol + col;
       postureIDs(i).postureID = postureID;
       postureIDs(i).workspaceCenter = workspaceCenter;      
    end
    [~,sortInd] = sort([postureIDs.postureID]);
    postureIDs = postureIDs(sortInd);
    
%% Keep only gridReachingData; get successful trials
    Data = gridReachData;
    clearvars gridReachData
    sucData = Data([Data.trialStatus]==1);

%% Get Reach Trajectories
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trajFields = {'marker'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'GridReaching'};
        trialInclStates(1).inclStates = {'Target Acquire'}; trialInclStates(1).inclOccurrence = {'first'};
%         trialInclStates(1).addTimeToBeginning = {0};
%         trialInclStates(1).addTimeToEnd = {0};
    trajStruct = getTrajStruct(sucData,condFields,trajFields,trialInclStates);    

%% Get Kinematic Data
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'},{'delayLength','kinData','delayLength'}};
    kinFields = {'reactTime','acqTime'};
    kinStruct = getKinStruct(sucData,condFields,kinFields);
    
%% Get Condition Lists
    kinData = [Data.kinData];
    targetList = unique([trajStruct.target]);
    postureList = unique([trajStruct.posture]);
    delayList = unique([kinStruct.delayLength]);
    holdTimeList = unique([kinData.targetHoldTime]);
    
%% Get P+T Stats
    PTStatsStruct = struct('posture',[],'target',[],'total',[],'suc',[],'sucPct',[],'longTHTotal',[],'longTHSuc',[],'longTHSucPct',[],'shortDelayTotal',[],'shortDelaySuc',[],'shortDelaySucPct',[],'longDelayTotal',[],'longDelaySuc',[],'longDelaySucPct',[],'medianAcqTime',[]);
    
    %Get raw data
    targetData = [Data.targetData];
    targetID = [targetData.targetID];
    conditionData = [Data.conditionData];
    postureID = [conditionData.postureID];
    trialStatus = double([Data.trialStatus]);
    targetHoldTime = [kinData.targetHoldTime];
    delayLength = [kinData.delayLength];
    centerHoldTime = [kinData.centerHoldTime];
    
    %Fill Struct
    structInd = 1;
    for posture = postureList
        for target = targetList
            PTStatsStruct(structInd).posture = posture;
            PTStatsStruct(structInd).target = target;
            %PT Total
            sucData = trialStatus([targetID==target]&[postureID==posture]);
            totalTrials = size(sucData,2);
            sucTrials = sum(sucData);
            sucPct = 100*sucTrials/totalTrials;
            PTStatsStruct(structInd).total = totalTrials;
            PTStatsStruct(structInd).suc = sucTrials;
            PTStatsStruct(structInd).sucPct = sucPct;
            %Long TH
            sucData = trialStatus([targetID==target]&[postureID==posture]&[targetHoldTime==1000]);
            totalTrials = size(sucData,2);
            sucTrials = sum(sucData);
            sucPct = 100*sucTrials/totalTrials;
            PTStatsStruct(structInd).longTHTotal = totalTrials;
            PTStatsStruct(structInd).longTHSuc = sucTrials;
            PTStatsStruct(structInd).longTHSucPct = sucPct;
            %Short Delay
            sucData = trialStatus([targetID==target]&[postureID==posture]&[delayLength<500]);
            totalTrials = size(sucData,2);
            sucTrials = sum(sucData);
            sucPct = 100*sucTrials/totalTrials;
            PTStatsStruct(structInd).shortDelayTotal = totalTrials;
            PTStatsStruct(structInd).shortDelaySuc = sucTrials;
            PTStatsStruct(structInd).shortDelaySucPct = sucPct;
            %Long Delay
            sucData = trialStatus([targetID==target]&[postureID==posture]&[delayLength>=500]);
            totalTrials = size(sucData,2);
            sucTrials = sum(sucData);
            sucPct = 100*sucTrials/totalTrials;
            PTStatsStruct(structInd).longDelayTotal = totalTrials;
            PTStatsStruct(structInd).longDelaySuc = sucTrials;
            PTStatsStruct(structInd).longDelaySucPct = sucPct;
            %Acq Time
            medianAcqTime = median([kinData([targetID==target]&[postureID==posture]).acqTime]);
            PTStatsStruct(structInd).medianAcqTime = medianAcqTime;
            %Update Struct Ind
            structInd = structInd +1;
        end
    end
    
    longCHSuc = struct('posture',[],'total',[],'suc',[],'sucPct',[]);
    structInd = 1;
    for posture = postureList
        sucData = trialStatus([postureID==posture]&[centerHoldTime==1000]);
        totalTrials = size(sucData,2);
        sucTrials = sum(sucData);
        sucPct = 100*sucTrials/totalTrials;
        longCHSuc(structInd).total = totalTrials;
        longCHSuc(structInd).suc = sucTrials;
        longCHSuc(structInd).sucPct = sucPct;
        structInd = structInd + 1;
    end
    
%% Add kinematic means, CI's 
    numCond = size(kinStruct,2);
    for i = 1:size(kinStruct,2)
       allReactTime = kinStruct(i).allReactTime;
       kinStruct(i).avgReactTime = mean(allReactTime);
       kinStruct(i).CIReactTime = std(allReactTime);
    end
    
%% Define color maps
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    cmap = customRainbow;
    hsvColor = rgb2hsv(cmap);
    hsvColor(:,2)=.3;
    cmap2 =hsv2rgb(hsvColor);
    
%% Plot Reach Trajectories and P+T Success Rates
    f = figure; f.Position = [10 10 1400 800];
    postureInd = 1;

    for posture = postureList
        subplot(3,3,postureInd)
        hold on
        for target = targetList
            tempData = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target);
            if ~isempty(tempData)
                numTrials = size(tempData.allMarker,2);
                for trial = 1:numTrials
                    traj = tempData.allMarker(trial).traj;
                    plot(traj(:,1),traj(:,2),'Color',cmap2(target,:),'LineWidth',0.5)
                    hold on
                end
                avgTraj = tempData.avgMarker.traj;
                plot(avgTraj(:,1),avgTraj(:,2),'Color',cmap(target,:),'LineWidth',3);
            end
            targetAngle = (target-1)*45;
            %Get success info
            PTStats = PTStatsStruct([PTStatsStruct.posture]==posture & [PTStatsStruct.target]==target);
            
            %
            textDist = 100;
            textX = textDist*cosd(targetAngle);
            textY = textDist*sind(targetAngle);
            
            if target ~= 7
                text(textX,textY,[num2str(round(PTStats.sucPct)),'% ','(',num2str(PTStats.suc),'/',num2str(PTStats.total),')',', MAT = ',num2str(round(PTStats.medianAcqTime))...
                    newline,'LTH ',num2str(round(PTStats.longTHSucPct)),'% ','(',num2str(PTStats.longTHSuc),'/',num2str(PTStats.longTHTotal),')',...
                    newline,'SD ',num2str(round(PTStats.shortDelaySucPct)),'% ','(',num2str(PTStats.shortDelaySuc),'/',num2str(PTStats.shortDelayTotal),')',...
                    ', LD ',num2str(round(PTStats.longDelaySucPct)),'% ','(',num2str(PTStats.longDelaySuc),'/',num2str(PTStats.longDelayTotal),')',...
                    ],'FontSize',8)
            else
                text(textX,textY,[num2str(round(PTStats.sucPct)),'% ','(',num2str(PTStats.suc),'/',num2str(PTStats.total),')',', MAT = ',num2str(round(PTStats.medianAcqTime))...
                    newline,'LTH ',num2str(round(PTStats.longTHSucPct)),'% ','(',num2str(PTStats.longTHSuc),'/',num2str(PTStats.longTHTotal),')',...
                    'SD ',num2str(round(PTStats.shortDelaySucPct)),'% ','(',num2str(PTStats.shortDelaySuc),'/',num2str(PTStats.shortDelayTotal),')',...
                    ', LD ',num2str(round(PTStats.longDelaySucPct)),'% ','(',num2str(PTStats.longDelaySuc),'/',num2str(PTStats.longDelayTotal),')',...
                    ],'FontSize',8)
            end
            
        end
        axis equal
        xlabel('x (mm)')
        ylabel('y (mm)')
        xlim([-100 100])
        ylim([-100 100])
        xticks([-100 100])
        yticks([-100 100])
        postureInd = postureInd + 1;
    end
    
    if saveFig
        saveas(gcf,[savePath,'reachTraj',dateStr,'.jpg'])
    end
    
%% Plot reaction time by delay length
    figure
    meanRxnTime = nan(1,size(delayList,2));
    delayInd = 1;
    for delay = delayList
       tempData = kinStruct([kinStruct.delayLength]==delay);
       allReactTime = [tempData.allReactTime];
       allReactTime(allReactTime > 1000 | allReactTime < 125) = [];
       for i = 1:size(allReactTime,2)
           jitter = 2*randn;
           plot(delay+jitter,allReactTime(i),'.','MarkerSize',10,'Color','k');
           hold on
       end
       meanRxnTime(delayInd) = mean(allReactTime);
       delayInd = delayInd + 1;
    end
%     for target = targetList
%         tempData = kinStruct([kinStruct.target]==target);
%         plot([tempData.delayLength],[tempData.avgReactTime],'Color',cmap(target,:))
%     end
    plot(delayList,meanRxnTime,'Color','k','LineWidth',3)
    xticks(delayList)
    xlabel('Delay (ms)')
    ylabel('Reaction Time (ms)')
    
    if saveFig
        saveas(gcf,[savePath,'reactTime',dateStr,'.jpg'])
    end
    
%% Get reach time distribution
    figure
    histogram([kinData.acqTime])
    xlabel('Reach Time (ms)')
    ylabel('Count')
    
    if saveFig
        saveas(gcf,[savePath,'reachTime',dateStr,'.jpg'])
    end
    
%% Plot workspace info
figure
hold on
for posture = postureList
   tempData = Data(postureID == posture);
   workspaceCenter = tempData(end).targetData.workspaceCenter;
   centerSize = tempData(end).targetData.centerSize;
   circle(workspaceCenter(1),workspaceCenter(2),centerSize,'k');
   for target = targetList
       tempData = Data(postureID == posture & targetID == target);
       targetLoc = tempData(end).targetData.targetLoc+workspaceCenter;
       targetSize = tempData(end).targetData.targetSize;
       circle(targetLoc(1),targetLoc(2),targetSize,'k');
   end
end
axis equal
xlabel('x (mm)')
ylabel('y (mm)')
    if saveFig
        saveas(gcf,[savePath,'workspace',dateStr,'.jpg'])
    end
% % Plot success rate by delay and hold 
%     figure 
%     subplot(2,4,4)
%     plot([delaySuc.condition],[delaySuc.pct])
%     xlabel('Delay (ms)')
%     ylabel('Success Rate (%)')
%     ylim([0 100])
%     
%     figure 
%     subplot(2,4,8)
%     plot([holdSuc.condition],[holdSuc.pct])
%     xlabel('Target Hold (ms)')
%     ylabel('Success Rate (%)')
%     ylim([0 100])
%     
% %% Title figure and Save
%     sgtitle(date)
%     if saveFig 
%         saveas(gcf,[savePath,'DCOdailyAnalysis',date,'.jpg'])
%     end