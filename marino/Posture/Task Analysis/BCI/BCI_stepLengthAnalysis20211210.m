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
                step1Times = [kinData.step1AcqTime];
                step1Times = step1Times(step1Times<900);
                step1Times = median(step1Times);
                jitter = 0.25*rand(1,length(step1Times));
                plot(step1Times,posture*ones(1,length(step1Times))+jitter,'.','MarkerSize',15,'Color',tcmap(target,:));
%                 text(rxnTimes,posture*ones(1,length(rxnTimes))+jitter,num2str(trialNum'));
            end
        end
    end
    xlabel('Median Step 1 Time (ms)')
    ylabel('Posture Label')
    
%% Plot reaction time histogram
    kinData = [Data.kinData];
    rxnTimes = [kinData.rxnTime];
    figure
    histogram(rxnTimes);
    xlabel('Reaction Time (ms)')
    ylabel('Count')
    
%% Plot movement time histogram
    kinData = [Data.kinData];
    movementTimes = [kinData.movementTime];
    figure
    histogram(movementTimes);
    xlabel('Movement Time (ms)')
    ylabel('Count')
    
%% Plot movement times by condition
    figure
    hold on
    for posture = 1:5
        for target = 1:8
            trialNum = [Data([conditionData.postureID]==posture & [targetData.target1ID]==target).trialNum];
            kinData = [Data([conditionData.postureID]==posture & [targetData.target1ID]==target).kinData];
            if ~isempty(kinData)
                movementTimes = [kinData.movementTime];
                jitter = 0.25*rand(1,length(movementTimes));
                plot(movementTimes,posture*ones(1,length(movementTimes))+jitter,'.','MarkerSize',7,'Color',tcmap(target,:));
%                 text(rxnTimes,posture*ones(1,length(rxnTimes))+jitter,num2str(trialNum'));
            end
        end
    end
    xlabel('Movement Time (ms)')
    ylabel('Posture Label')
    
    
%% Plot median step 1 times by condition
    figure
    hold on
    for posture = 1:5
        for target = 1:8
            trialNum = [Data([conditionData.postureID]==posture & [targetData.target1ID]==target).trialNum];
            kinData = [Data([conditionData.postureID]==posture & [targetData.target1ID]==target).kinData];
            if ~isempty(kinData)
                movementTimes = [kinData.movementTime];
                movementTimes = movementTimes(movementTimes<900);
                movementTimes = median(movementTimes);
                jitter = 0.25*rand(1,length(movementTimes));
                plot(movementTimes,posture*ones(1,length(movementTimes))+jitter,'.','MarkerSize',15,'Color',tcmap(target,:));
%                 text(rxnTimes,posture*ones(1,length(rxnTimes))+jitter,num2str(trialNum'));
            end
        end
    end
    xlabel('Median Movement Time (ms)')
    ylabel('Posture Label')