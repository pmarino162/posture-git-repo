   clear; clc; clf; close all;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = customRainbow;
    
%% Load data
    [Data] = loadEarlData20200116_20211210();
    
%% Get trajStruct
    trajFields = {'forceCursor'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;
    trialInclStates(1).trialName = {'IsometricForce_1D'};
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).inclStates = {{'state','Center Hold','last',0},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    
%% Plot targets and trajectories
    targetData = [Data.targetData];
    targetID = [targetData.targetID];
    targetList = unique(targetID);
    %Force cursor fig
    figure; hold on;
    for posture = 1
        for target = targetList
            %Plot target
                trialInd = find([targetID==target],1);
                targetLoc = Data(trialInd).targetData.targetLoc;
                targetSize = Data(trialInd).targetData.targetSize;
                plot(targetLoc(1,1),targetLoc(1,2),'.','MarkerSize',20,'Color',tcmap(target,:));
            %Plot cursor traj
                tempTrajStruct = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture);
%                 numTraj = size(tempTrajStruct.allForceCursor,2);
%                 for i = 1:numTraj
%                    traj = tempTrajStruct.allForceCursor(i).traj; 
%                    timestamps = tempTrajStruct.allForceCursor(i).timestamps; 
%                    plot(traj(:,1),traj(:,2),'Color',tcmap(target,:))
%                 end
                traj = tempTrajStruct.avgForceCursor.traj;
                plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
        end
    end
    xlabel('x (mm)')
    ylabel('y (mm)')
    
    
    %Force cursor v time
    figure;
    for posture = 1
        for target = targetList
            %Plot cursor traj
                tempTrajStruct = trajStruct([trajStruct.target]==target & [trajStruct.posture]==posture);
                traj = tempTrajStruct.avgForceCursor.traj;
                time = tempTrajStruct.avgForceCursor.timestamps;
                subplot(2,1,1)
                 hold on;
                    plot(time,traj(:,1),'Color',tcmap(target,:))
                    xlabel('t')
                    ylabel('x (mm)')
                subplot(2,1,2)
                 hold on;
                    plot(time,traj(:,2),'Color',tcmap(target,:))
                    ylabel('y (mm)')
                    xlabel('t')
        end
    end
   

    %Force fig
    
%% Try rxn time algorithm
    %Get trialData
    trial = 1;
    trialData = Data(trial);
    trialName = trialData.trialName;
    
    %Set up trial include states
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {trialName};
    trialInclStates(1).addTimeToBeginning = {0};
    trialInclStates(1).addTimeToEnd = {0};
    %Get state data
    stateNames = trialData.stateData.stateNames;
    stateTransitions = trialData.stateData.stateTransitions;
    targetID = max(find([cellfun(@(x) strcmpi(x,'Target'),stateNames)]==1));
    targetHoldID = max(find([cellfun(@(x) strcmpi(x,'Target Hold'),stateNames)]==1));
    successID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
    %Get state transitions and step acquisition times
    if ismember(targetID,stateTransitions(1,:))
        targetTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetID))));
    end
    if ismember(targetHoldID,stateTransitions(1,:))
        targetHoldTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetHoldID))));
        targetAcqTime = targetHoldTime - targetTime;
        kinData.targetAcqTime = targetAcqTime;
    else
        kinData = [];
    end
    %Get movement onset, reaction, and movement times if trial made it to
    %target hold
    if ismember(targetHoldID,stateTransitions(1,:))
        %Get baseline force distribution (before go cue)
        trialInclStates(1).inclStates = {{'state','Center Hold','last',0},{'state','Target','first',0}};
        [centerForceTraj,~] = getStatesTraj20211210(trialData,trialInclStates,'force',1,'timeRelToTrialStart',true);
        meanYForce = mean(centerForceTraj(:,2));
        stdYForce = std(centerForceTraj(:,2));
        %Get Y force trajectory after go cue
        trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
        [targetForceTraj,timestamps] = getStatesTraj20211210(trialData,trialInclStates,'force',1,'timeRelToTrialStart',true);
        YTargetForceTraj = targetForceTraj(:,2);
        %Walk backwards from target hold time to nearest time when force
        %went 20% towards end force above baseline mean in correct direction. Interpolate to improve estimate
        threshold = 0.2*(YTargetForceTraj(end)-meanYForce) + meanYForce;
        i = length(YTargetForceTraj);
        if threshold > 0
            while YTargetForceTraj(i) > threshold
               i = i-1; 
            end
            moveOnsetTime = interp1(YTargetForceTraj(i:i+1),timestamps(i:i+1),threshold);
        elseif threshold < 0 
            while YTargetForceTraj(i) < threshold
               i = i-1; 
            end
            moveOnsetTime = interp1(YTargetForceTraj(i:i+1),timestamps(i:i+1),threshold);
        end
        kinData.moveOnsetTime = moveOnsetTime;
        kinData.rxnTime = moveOnsetTime - targetTime;
    end
    
%% Plot rxnTime Histogram
    kinData = [Data.kinData];
    rxnTime = [kinData.rxnTime];
    histogram(rxnTime)
    
    targetData = [Data.targetData];
    targetID = [targetData.targetID];
    
    t3rxnTimes = rxnTime(targetID==3);
    t7rxnTimes = rxnTime(targetID==7);
    figure; hold on;
    histogram(t3rxnTimes);
    histogram(t7rxnTimes);