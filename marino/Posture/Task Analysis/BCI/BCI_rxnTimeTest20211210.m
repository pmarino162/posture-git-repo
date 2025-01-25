clear; clc; clf; close all;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = customRainbow;
    
%% Set up save 
    saveFig = true;
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20211229 - adam meeting - controlling for kin and are traj parallel';

%% Load Data 
   [Data] = loadEarlData20200317_20211210;
   conditionData = [Data.conditionData];
   targetData = [Data.targetData];
   
%% Test rxnTime identification procedure  
   trial = 1;
    trialData = Data(trial);
    trialName = 'GridTask_CO_Across_BC_ForceBar';
    %Set up trial include states
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {trialName};
    trialInclStates(1).addTimeToBeginning = {0};
    trialInclStates(1).addTimeToEnd = {0};
    trialInclStates(1).inclStates = {'Step 1'};
    trialInclStates(1).inclOccurrence = {'first'};


   centerLoc = Data(trial).targetData.centerLoc;
   target1Loc = Data(trial).targetData.target1Loc;
   target1Vec = target1Loc(1,1:2)-centerLoc(1,1:2);
   target1Dist = norm(target1Vec);
   [cursorTraj,timestamps,~] = getStatesTraj(trialData,trialInclStates,'decoderCursorTraj','timeRelToTrialStart',true);
   target1Proj = dot(cursorTraj,repmat(target1Vec,length(timestamps),1),2)./target1Dist;
   threshInd = min(find(target1Proj > 0.5*target1Dist));
   moveOnsetTime = interp1(target1Proj(threshInd-1:threshInd),timestamps(threshInd-1:threshInd),0.5*target1Dist);
   
   
%% Plot step 1 times by condition
    figure
    hold on
    for posture = 1:5
        for target = 1:8
            trialNum = [Data([conditionData.postureID]==posture & [targetData.target1ID]==target).trialNum];
            kinData = [Data([conditionData.postureID]==posture & [targetData.target1ID]==target).kinData];
            if ~isempty(kinData)
                step1Times = [kinData.step1AcqTime];
                jitter = 0.25*rand(1,length(step1Times));
                plot(step1Times,posture*ones(1,length(step1Times))+jitter,'.','MarkerSize',7,'Color',tcmap(target,:));
%                 text(rxnTimes,posture*ones(1,length(rxnTimes))+jitter,num2str(trialNum'));
            end
        end
    end
    xlabel('Step 1 Time (ms)')
    ylabel('Posture Label')
   
%% Plot median step 1 times by condition
    figure
    hold on
    for posture = 1:5
        for target = 1:8
            trialNum = [Data([conditionData.postureID]==posture & [targetData.target1ID]==target).trialNum];
            kinData = [Data([conditionData.postureID]==posture & [targetData.target1ID]==target).kinData];
            if ~isempty(kinData)
                step1Times = median([kinData.step1AcqTime]);
                jitter = 0.25*rand(1,length(step1Times));
                plot(step1Times,posture*ones(1,length(step1Times))+jitter,'.','MarkerSize',7,'Color',tcmap(target,:));
%                 text(rxnTimes,posture*ones(1,length(rxnTimes))+jitter,num2str(trialNum'));
            end
        end
    end
    xlabel('Reaction Time (ms)')
    ylabel('Posture Label')