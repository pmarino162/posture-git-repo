function [Data] = getDataStructChaseLabDissociation20220419(rawData,varargin)

% This function processes raw experimental data so that it's in a
% user-friendly format

%% Variable Arguments
    getSpikes = false;
        getSorts = false;
        exclCh = [];
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = true;
        centerMarker = true; %Subtracts workspace center from marker pose if true
    inclStateTable = false;
    getKin = false;
    exclTrials = [];
    trialName = '';
    assignopts(who,varargin);
    
%% Preallocate Data Struct
    %stateData
    stateData = struct('stateNames',{},'stateTransitions',[]);
    %Data
    Data = struct('trialNum',[],'trialName','','trialStatus',[],'targetData',[],'conditionData',[],'stateData',stateData,'cursor',[]);
    
    %Marker
    if getMarker == true
        marker = struct('position',[],'velocity',[],'acceleration',[]);
        Data.marker = marker;
    end
    %Spikes
    if getSpikes == true
        %Data.spikes = struct('channel',[],'sort',[],'timestamps',[]);
    end

    taskInfo = rawData.taskInfo;
    rawData = rawData.trialData;
    numTrials = size(rawData,1);
    Data = repmat(Data,1,numTrials);

%% Fill Data Struct
for trial = 1:numTrials
    %Get Trial Number
        Data(trial).trialNum = rawData(trial).trial;
    %Get Trial Name
        Data(trial).trialName = trialName;
    %Get time
        Data(trial).time = rawData(trial).time;
    %Get Trial Status
        Data(trial).trialStatus = rawData(trial).trialStatusLabel;
    %Get Target Data
        Data(trial).targetData.centerLoc = rawData(trial).z_extraInfo.startTargetPosition.*1000; %Convert to mm
        Data(trial).targetData.centerSize = taskInfo.targetRadius;
        Data(trial).targetData.targetLoc = rawData(trial).z_extraInfo.endTargetPosition.*1000; %Convert to mm
        Data(trial).targetData.targetSize = taskInfo.targetRadius;
        Data(trial).targetData.targetID = rawData(trial).directionLabel;
        Data(trial).targetData.cursorRadius = taskInfo.cursorRadius;
        Data(trial).targetData.workspaceCenter = Data(trial).targetData.centerLoc;
    %Get Condition Data
        Data(trial).conditionData.postureID = rawData(trial).proprioLabel;
        Data(trial).conditionData.visualID = rawData(trial).visualLabel;
    %Get State Data (and convert transitions to ms)
        stateNames = {'Trial Start','Center Hold','Reach','Target Hold','Success with Reward','Trial End'};
        stateTransitions = zeros(2,6);
            stateTransitions(1,:) = 1:6;
            stateTransitions(2,1) = rawData(trial).time(1);
            stateTransitions(2,2) = rawData(trial).t1_centerHoldStart;
            stateTransitions(2,3) = rawData(trial).t2_goCueTime;
            stateTransitions(2,4) = rawData(trial).t6_postReachStateTime;
            stateTransitions(2,5) = rawData(trial).t7_trialEndTime;
            stateTransitions(2,6) = rawData(trial).time(end);
        Data(trial).stateData(1).stateNames = stateNames;
        Data(trial).stateData(1).stateTransitions = stateTransitions;
    %Get Spike Data
        if getSpikes == true           
            Data(trial).spikes = rawData(trial).spikeMatrix';
        end
    %Get Marker Data (Phasespace)
        if getMarker == true
            Data(trial).marker.position = rawData(trial).kinematics_updated.marker.position;
            Data(trial).marker.velocity = rawData(trial).kinematics_updated.marker.velocity;
            Data(trial).marker.acceleration = rawData(trial).kinematics_updated.marker.acceleration;
        end
    %Get Cursor Data
        Data(trial).cursor = rawData(trial).kinematics_updated.cursor;
    %Get Kin Data - REPLACED WITH FXN AT END OF
    %preprocessAndSaveData20220419
        %Data(trial).kinData.rxnTime = rawData(trial).t3_reactionTime;
        %Data(trial).kinData.peakSpeedTime = rawData(trial).t4_peakSpeedTime;
        %Data(trial).kinData.reachEndTime = rawData(trial).t5_reachEndTime;
end
     
end