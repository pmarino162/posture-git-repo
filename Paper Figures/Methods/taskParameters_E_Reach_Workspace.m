clear; clc;

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;

%% Load data, getTrajStruct
    dataset = 'E20210706';

    [Data,zScoreParams] = loadData(dataset);
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    trajFields = {'zSmoothFR','marker'};
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 

%% Look at individual trials
    postureID = [Data.conditionData];
    postureID = [postureID.postureID];
    targetID = [Data.targetData];
    targetID = [targetID.targetID];
    
    %Posture 1 Target 1
    posture = 1;
    target = 1;
    trial = min(find(postureID==posture & targetID==target));
    centerLoc = Data(trial).targetData.centerLoc
    targetLoc = Data(trial).targetData.targetLoc
    workspaceCenter = Data(trial).targetData.workspaceCenter
    
    %Posture 7 Target 1
    posture = 7;
    target = 1;
    trial = min(find(postureID==posture & targetID==target))
    centerLoc = Data(trial).targetData.centerLoc
    targetLoc = Data(trial).targetData.targetLoc
    workspaceCenter = Data(trial).targetData.workspaceCenter
    
%% Plot Positions using trajStruct 
    figure
    for posture = postureList
        for target = targetList
            if ~isempty(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target))
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj;
                if posture == 1
                    plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
                else
                    plot(traj(:,1),traj(:,2),'--','Color',tcmap(target,:));
                end
                hold on
            end
        end
    end