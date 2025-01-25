clear; clc; clf; close all;

%% Setup Save
    savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Task Analysis\Grid Reaching\';
    saveFig = false;    
    
%% Load Data 
    %Earl - 7/06/2021
    date = '20210706';
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
    % Get Data Struct, 
    binWidth = 25;
    exclCh =  [44 87 88 77 78 71 67 69 118];
    exclCh =  [];
    getSorts = false;
    Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'binWidth',binWidth,'exclZero',false,'getMarker',true,'centerMarker',false,'getKin',true);
    % Clean Data
    [Data] = cleanData20210706(Data);
    %Keep only successful trials
    Data = Data([Data.trialStatus]==1);
    % Label Postures
    [Data,postureIDs] = labelPostures20210706(Data);
    %Keep select postures
    allTrialPostures = [Data.conditionData];
    allTrialPostures = [allTrialPostures.postureID];
    keepPostures = [1,2,4,5,8,12,14];
    Data = Data(ismember(allTrialPostures,keepPostures));
    % Keep trials appropriate for reach and for delay
    kinData = [Data.kinData];
    delayLength = [kinData.delayLength];
    delayData = Data(delayLength >= 500);
    
%% Get Reach Trajectories
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    if getSorts
        trajFields = {'allSortSmoothedFR','marker'};
    else
        trajFields = {'allChannelSmoothedFR','marker'};
    end
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'GridReaching'};
        trialInclStates(1).inclOccurrence = {'first'};
    %Delay
    trialInclStates(1).inclStates = {'Delay'}; 
    trialInclStates(1).addTimeToBeginning = {-100};
    delayTrajStruct = getTrajStruct(delayData,condFields,trajFields,trialInclStates,'matchConditions',true);    
    %Reach 
    trialInclStates(1).inclStates = {'Target Acquire'}; 
    trialInclStates(1).addTimeToBeginning = {-100};
    reachTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,'matchConditions',true);

    % Keep only trials with length w/in 2 stdDevs of mean length
%     [Data] = excludeLengths(Data,trialInclStates);
    
%% Define color maps
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
%     cmap = customRainbow;
%     hsvColor = rgb2hsv(cmap);
%     hsvColor(:,2)=.3;
%     cmap2 =hsv2rgb(hsvColor);    

%% Plot PSTH
    figure
    channel = 119;
    cmap = customRainbow;
    reachTimeOffset = 1100;
    for target = 1:8
        posture = 1;
        %Delay
        tempData = delayTrajStruct([delayTrajStruct.target]==target & [delayTrajStruct.posture]==posture);
        
        traj = tempData.avgAllChannelSmoothedFR.traj(:,channel);
        timestamps = tempData.avgAllChannelSmoothedFR.timestamps;
        plot(timestamps,traj,'Color',cmap(target,:),'LineWidth',1)
        
        hold on
        %Reach
        tempData = reachTrajStruct([reachTrajStruct.target]==target & [reachTrajStruct.posture]==posture);
        
        traj = tempData.avgAllChannelSmoothedFR.traj(:,channel);
        timestamps = tempData.avgAllChannelSmoothedFR.timestamps+reachTimeOffset;
        plot(timestamps,traj,'Color',cmap(target,:),'LineWidth',1)
    end
    xlabel('time (ms)')
    ylabel('FR (Hz)')
    xlim([-250 reachTimeOffset+1000])
    xticks([-250 0 500 reachTimeOffset reachTimeOffset+500])
    xticklabels({'-250','TO','500','MO','500'})
%     xticklabels({'TO',})
%     figure
%     cmap = customSummer;
%     for posture = [1,2,4,5]
%         target = 1;
%         tempData = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture);
%         traj = tempData.avgAllChannelSmoothedFR.traj(:,channel);
%         timestamps = tempData.avgAllChannelSmoothedFR.timestamps;
%         plot(timestamps,traj,'Color',cmap(posture,:),'LineWidth',1)
%         hold on
%     end
%     xlabel('time (ms)')
%     ylabel('FR (Hz)')