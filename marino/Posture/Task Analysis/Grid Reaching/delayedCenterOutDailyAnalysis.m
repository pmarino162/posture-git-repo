clear; clc; clf; close all;

%% Load Data
    date = '20210622';
    load('D:\Animals\Earl\2021\06\20210622\Earl20210622_Delayed Center Out_translated.mat')
    Data = getDataStruct(Data,'getMarker',true,'getKin',true);
    sucData = Data([Data.trialStatus]==1);

%% Get Reach Trajectories
    condFields = {{'target','targetData','targetID'}};
    trajFields = {'marker'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'Delayed Center Out 20210621'};
        trialInclStates(1).inclStates = {'Target Acquire'}; trialInclStates(1).inclOccurrence = {'first'};
        trialInclStates(1).addTimeToBeginning = {0};
        trialInclStates(1).addTimeToEnd = {0};
    trajStruct = getTrajStruct(sucData,condFields,trajFields,trialInclStates);    

%% Get Kinematic Data
    condFields = {{'target','targetData','targetID'},{'delayLength','kinData','delayLength'}};
    kinFields = {'reactTime','acqTime'};
    kinStruct = getKinStruct(sucData,condFields,kinFields);
    
%% Get Condition Info
    targetList = unique([trajStruct.target]);
    delayList = unique([kinStruct.delayLength]);
    kinData = [Data.kinData];
    holdTimeList = unique([kinData.targetHoldTime]);
    
%% Get Success Rates
    sucStruct = struct('condition',[],'total',[],'suc',[],'pct',[]);
    %by target
    numTargets = size(targetList,2);
    targSuc = repmat(sucStruct,1,numTargets);
    targetData = [Data.targetData];
    targetID = [targetData.targetID];
    structInd = 1;
    for target = targetList
        tempData = Data(targetID == target);
        total = size(tempData,2);
        suc = sum(double([tempData.trialStatus]));
        targSuc(structInd).condition = target;
        targSuc(structInd).total = total;
        targSuc(structInd).suc = suc;
        targSuc(structInd).pct = (suc/total)*100;
        structInd = structInd + 1;
    end
    %by delay
    numDelays = size(delayList,2);
    delaySuc = repmat(sucStruct,1,numDelays);
    kinData = [Data.kinData];
    delayLength = [kinData.delayLength];
    structInd = 1;
    for delay = delayList
        tempData = Data(delayLength == delay);
        total = size(tempData,2);
        suc = sum(double([tempData.trialStatus]));
        delaySuc(structInd).condition = delay;
        delaySuc(structInd).total = total;
        delaySuc(structInd).suc = suc;
        delaySuc(structInd).pct = (suc/total)*100;
        structInd = structInd + 1;
    end
    %by holdTime
    numHoldTimes = size(holdTimeList,2);
    holdTimeSuc = repmat(sucStruct,1,numHoldTimes);
    kinData = [Data.kinData];
    targetHoldTime = [kinData.targetHoldTime];
    structInd = 1;
    for holdTime = holdTimeList
        tempData = Data(holdTime == targetHoldTime);
        total = size(tempData,2);
        suc = sum(double([tempData.trialStatus]));
        holdSuc(structInd).condition = holdTime;
        holdSuc(structInd).total = total;
        holdSuc(structInd).suc = suc;
        holdSuc(structInd).pct = (suc/total)*100;
        structInd = structInd + 1;
    end
    
%% Add kinematic means, CI's 
    numCond = size(kinStruct,2);
    for i = 1:size(kinStruct,2)
       allReactTime = kinStruct(i).allReactTime;
       kinStruct(i).avgReactTime = mean(allReactTime);
       kinStruct(i).CIReactTime = std(allReactTime);
    end
    
%% Define figure and color maps
    f = figure; f.Position = [10 10 1400 800];
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    cmap = customRainbow;
    hsvColor = rgb2hsv(cmap);
    hsvColor(:,2)=.3;
    cmap2 =hsv2rgb(hsvColor);
    
%% Plot Reach Trajectories
    subplot(2,4,[1,2,5,6])
    hold on
    for target = targetList
        tempData = trajStruct([trajStruct.target]==target);
        numTrials = size(tempData.allMarker,2);
        for trial = 1:numTrials
            traj = tempData.allMarker(trial).traj;
            plot(traj(:,1),traj(:,2),'Color',cmap2(target,:),'LineWidth',0.5)
            hold on
        end
        avgTraj = tempData.avgMarker.traj;
        plot(avgTraj(:,1),avgTraj(:,2),'Color',cmap(target,:),'LineWidth',3);
        targetAngle = (target-1)*45;
        textDist = 90;
        textX = textDist*cosd(targetAngle);
        textY = textDist*sind(targetAngle);
        text(textX,textY,[num2str(round(targSuc(target).pct)),'%'],'FontSize',14)
    end
    axis equal
    xlabel('x (mm)')
    ylabel('y (mm)')
    xlim([-100 100])
    ylim([-100 100])
    xticks([-100 100])
    yticks([-100 100])
    
%% Plot reaction time by delay length
%     figure
subplot(2,4,3)
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
    for target = targetList
        tempData = kinStruct([kinStruct.target]==target);
        plot([tempData.delayLength],[tempData.avgReactTime],'Color',cmap(target,:))
    end
    plot(delayList,meanRxnTime,'Color','k','LineWidth',2)
    xlabel('Delay (ms)')
    ylabel('Reaction Time (ms)')
    
%% Get reach time distribution
%     figure
    subplot(2,4,7)
    histogram([kinData.acqTime])
    xlabel('Reach Time (ms)')
    ylabel('Count')
    
%% Plot success rate by delay and hold 
%     figure 
    subplot(2,4,4)
    plot([delaySuc.condition],[delaySuc.pct])
    xlabel('Delay (ms)')
    ylabel('Success Rate (%)')
    ylim([0 100])
    
%     figure 
    subplot(2,4,8)
    plot([holdSuc.condition],[holdSuc.pct])
    xlabel('Target Hold (ms)')
    ylabel('Success Rate (%)')
    ylim([0 100])
    
    sgtitle(date)