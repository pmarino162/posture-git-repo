clear; clc; clf; close all;

%% Choose Task
%     task = 'BCI';
%     task = 'Iso';
    task = 'Reach';
    
%% Setup save dir   
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20210930 - first model\';

%% Load Data 
    binWidth = 25;
    switch task
        case 'BCI'    
            [Data,~,~,~,~,...
             ~,~,~,~,...
             ~,~] = loadEarlData20200317(binWidth);
        case 'Iso'
            [Data] = loadEarlData20200116(binWidth);
        case 'Reach'
            load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
            exclCh =  [44 87 88 77 78 71 67 69 118];
            getSorts = false;
            Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'binWidth',binWidth,'exclZero',false,'getMarker',true,'getKin',true);
            [Data] = cleanData20210706(Data);
            Data = Data([Data.trialStatus]==1);
            [Data,postureIDs] = labelPostures20210706(Data);
            allTrialPostures = [Data.conditionData];
            allTrialPostures = [allTrialPostures.postureID];
            keepPostures = [1,2,4,5];
%             keepPostures = [1,2,4,5,8,12,14];
            Data = Data(ismember(allTrialPostures,keepPostures));
    end

%% Get Traj Struct   
    trajFields = {'allChannelSmoothedFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    switch task
        case 'BCI'
            condFields = {{'posture','conditionData','postureID'},{'target','targetData','target1ID'}};
            trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
            trialInclStates(1).inclStates = {'Step 1'};
            trialInclStates(1).inclOccurrence = {'last'};
                trialInclStates(1).addTimeToBeginning = {0};
                trialInclStates(1).addTimeToEnd = {0};   
        case 'Iso'
            condFields = {{'posture','conditionData','postureID'},{'target','targetData','targetID'}};
            trialInclStates(1).trialName = {'IsometricForce_1D'};
            trialInclStates(1).inclStates = {'Target'};
            trialInclStates(1).inclOccurrence = {'last'};
                trialInclStates(1).addTimeToBeginning = {0};
                trialInclStates(1).addTimeToEnd = {0}; 
        case 'Reach'
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).trialName = {'GridReaching'};
            trialInclStates(1).inclStates = {'Target Acquire'};
            trialInclStates(1).inclOccurrence = {'first'};
                trialInclStates(1).addTimeToBeginning = {0};
                trialInclStates(1).addTimeToEnd = {0}; 
    end
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,'matchConditions',true);

%% Trim average trajectories to length of shortest one
    avgTrajLengths = nan(1,size(trajStruct,2));
    for i = 1:size(trajStruct,2)
        avgTrajLengths(i) = size(trajStruct(i).avgAllChannelSmoothedFR.traj,1);
    end
    minAvgTrajLength = min(avgTrajLengths);
    for i = 1:size(trajStruct,2)
        trajStruct(i).avgAllChannelSmoothedFR.traj = trajStruct(i).avgAllChannelSmoothedFR.traj(1:minAvgTrajLength,:);
        trajStruct(i).avgAllChannelSmoothedFR.timestamps = trajStruct(i).avgAllChannelSmoothedFR.timestamps(1,1:minAvgTrajLength);
    end
    
%% Get Postures, Targets, and Channels
    postureList = unique([trajStruct.posture]);
    targetList = unique([trajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgAllChannelSmoothedFR.traj,2);

%% dPCA
    combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
    margNames = {'Posture', 'Target', 'Condition-independent', 'P/T Interaction'};
    margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

    N = 119;
    P = numPostures;
    T = numTargets;
    timePts = size(trajStruct(1).avgAllChannelSmoothedFR.traj,1);
    
    firingRatesAverage = zeros(N,P,T,timePts);
    postureInd = 1;
    targetInd = 1;
    for posture = postureList
        for target = targetList
            firingRatesAverage(:,postureInd,target,:) = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgAllChannelSmoothedFR.traj';
        end
        postureInd = postureInd + 1;
    end
    
    time = trajStruct(1).avgAllChannelSmoothedFR.timestamps;
    timeEvents = [];
    tic
[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams);
toc

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);
