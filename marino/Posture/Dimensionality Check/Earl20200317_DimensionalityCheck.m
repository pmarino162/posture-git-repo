clear;
clc;
clf;
close all;

%% Load Data
    %3/17/2020
    date = '20200317';
    exclCh =  [3 12 16 20 31 73 85 92 94];
    binWidth = 1;
    [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
     N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200317(binWidth);
 
%% Collect States of Interest into Traj Field
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    numTrials = size(Data,2);
    for trial = 1:numTrials
        %Step 1
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
        [step1FR,step1FRTimestamps] = getStatesTraj(Data(trial),trialInclStates,'allChannelSmoothedFR'); 
        Data(trial).Traj.step1FR = step1FR; Data(trial).Traj.step1FRTimestamps = step1FRTimestamps;
        [step1Marker,step1MarkerTimestamps] = getStatesTraj(Data(trial),trialInclStates,'marker'); 
        Data(trial).Traj.step1Marker = step1Marker; Data(trial).Traj.step1MarkerTimestamps = step1MarkerTimestamps;
    end
    
%% Keep only trials with length w/in 2 stdDevs of mean length
    %Determine Lengths to keep
    numTrials = size(Data,2);
    step1Lengths = zeros(1,numTrials);
    for trial = 1:numTrials
        step1Timestamps = Data(trial).Traj.step1FRTimestamps;
        step1Lengths(trial) = size(step1Timestamps,2);
    end
    meanStep1Length = mean(step1Lengths);
    stdStep1Length = std(step1Lengths);
    minKeepLength = meanStep1Length - 2*stdStep1Length;
    if minKeepLength < 0
        minKeepLength = 0;
    end
    maxKeepLength = meanStep1Length + 2*stdStep1Length;
    %Plot distribution and keep lengths
    figure
        histogram(step1Lengths)
        hold on
        ca = gca;
        maxY = ca.YLim(2);
        line([minKeepLength,minKeepLength],[0 maxY],'Color','red')
        line([maxKeepLength,maxKeepLength],[0 maxY],'Color','red')
        xlabel('# Bins')
        ylabel('Count');
        title('Step 1 Lengths')
    %Remove trials
    keepMask = step1Lengths > minKeepLength & step1Lengths < maxKeepLength;
    Data = Data(keepMask);
    
%% Create trajStruct
    condFields = {{'posture','postureID'},{'target','targetData','target1ID'}};
    trajFields = {'step1FR','step1Marker'};
    trajStruct = getTrajStruct(Data,condFields,trajFields);
% 
%     trajStruct = struct('Condition','','PostureID','','TargetID',[],'allStep1FRTraj',[],'avgStep1FRTraj',[]);
%     structInd = 1;
%     conditionList = unique({Data.Condition});
%     for condition = conditionList
%        trajStruct(structInd).Condition = condition{1,1};
%        %Get condition data, preallocate
%        tempData = Data(strcmpi({Data.Condition},condition));
%        trajStruct(structInd).Posture = tempData(1).Posture;
%        trajStruct(structInd).TargetID = tempData(1).targetData.targetID;
%        numTrials = size(tempData,2); 
%        step1FRTrajLengths = zeros(1,numTrials);
%        %Get traj lengths
%        for trial = 1:numTrials
%            step1FRTraj = tempData(trial).Traj.step1FRTraj;
%            step1FRLength = size(step1FRTraj,1);
%            step1FRTrajLengths(trial) = step1FRLength;
%        end
%        %Get max and mean traj lengths, preallocate
%            step1FRMaxLength = max(step1FRTrajLengths);
%            step1FRMeanLength = round(mean(step1FRTrajLengths));
%            allstep1FRTraj = NaN(step1FRMaxLength,numDims,numTrials);
%        %Store all neural traj
%        for trial = 1:numTrials
%            step1FRTraj = tempData(trial).Traj.step1FRTraj;
%            step1FRLength = size(step1FRTraj,1);
%            allstep1FRTraj(1:step1FRLength,:,trial) = step1FRTraj;
%            trajStruct(structInd).allStep1FRTraj(trial).Traj = step1FRTraj;
%            trajStruct(structInd).allStep1FRTraj(trial).Timestamps = tempData(trial).Traj.step1FRTrajTimestamps;
%        end
%        %Take condition averages, store
%            step1FRConditionAvgTraj =  nanmean(allStep1FRTraj,3);
%            step1FRConditionAvgTraj = step1FRConditionAvgTraj(1:step1FRMeanLength,:);
%            trajStruct(structInd).avgStep1FRTraj = step1FRConditionAvgTraj;
%            trajStruct(structInd).step1FRTrajLengths = step1FRTrajLengths;
%        %Update Struct Ind
%        structInd = structInd + 1;
%     end
%     
%     clearvars conditionList tempData Step1FRTrajLengths  Step1FRTraj ...
%         Step1FRLength   Step1FRMaxLength Step1FRMeanLength...
%         allStep1FRTraj     avgStep1FRTraj...
%         Step1FRConditionAvgTraj 
%% PCA

%% FA

%% GPFA
