clear; clc;

%% Setup colormap    
    load('C:\Users\pmari\OneDrive\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;

%% Load data, getTrajStruct
    %Rocky
%     dataset = 'R20200221';
%     rawData = load('D:\Animals\Rocky\2020\02\20200221\Rocky_20200221_preprocessedData_20220512_165624.mat');
%     
%     %Nigel
    dataset = 'N20190226';
    rawData = load('D:\Animals\Nigel\2019\02\20190226\Nigel_20190226_preprocessedData_20220506_080540.mat');

    [Data,zScoreParams] = loadData(dataset);
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    trajFields = {'zSmoothFR','marker'};
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 
     

%% Get workspace centers and target distances
    %Kept only trials with visualID = 2
    
    trialData = [rawData.trialData];
    visualLabel = [trialData.visualLabel];
    proprioLabel = [trialData.proprioLabel];
    directionLabel = [trialData.directionLabel];
    
    workspaceCenter1 = trialData(min(find(visualLabel==2 & proprioLabel==1))).z_extraInfo.cursorTransform(1:2,4)*1000;
    workspaceCenter2 = trialData(min(find(visualLabel==2 & proprioLabel==2))).z_extraInfo.cursorTransform(1:2,4)*1000;
    
    xRange = abs(workspaceCenter1(1)-workspaceCenter2(1));
    yRange = abs(workspaceCenter1(2)-workspaceCenter2(2));
    xMin = min(workspaceCenter1(1),workspaceCenter2(1));
    yMin = min(workspaceCenter1(2),workspaceCenter2(2));
    workspaceCenter = [xMin+xRange/2;yMin+yRange/2];
    workspaceCenter1 = workspaceCenter1 - workspaceCenter;
    workspaceCenter2 = workspaceCenter2 - workspaceCenter;
    workspaceDist = norm(workspaceCenter1-workspaceCenter2);
    
    startTarget = trialData(1).z_extraInfo.startTargetPosition;
    endTarget = trialData(1).z_extraInfo.endTargetPosition;
    targetDist = norm(endTarget-startTarget)*1000;
    
    cursorRadius = rawData.taskInfo.cursorRadius;
    targetRadius = rawData.taskInfo.targetRadius;

    
% %% Print all raw data target locations
% 
%     rawDataTargetLocations = zeros(8,2,2);
%     rawDataCenterLocations= zeros(8,2,2);
%     
%     for posture = 1:2
%         posture
%         postureTrialData = trialData(find(visualLabel==2 & proprioLabel==posture));
%         extraInfo = [postureTrialData.z_extraInfo];
%         startTargetPosition = unique(vertcat(extraInfo.startTargetPosition),'rows')
%         endTargetPosition = unique(vertcat(extraInfo.endTargetPosition),'rows')
%     end

%% Print target location in marker coordinates
    for posture = 1:2
        posture
        cursorTransform = trialData(min(find(visualLabel==2 & proprioLabel==posture))).z_extraInfo.cursorTransform(1:2,4)*1000
        for target = 1:8
            tempData = trialData(min(find(visualLabel==2 & proprioLabel==posture & directionLabel==target))).z_extraInfo.endTargetPosition(1:2).*1000;
            tempData - cursorTransform'
        end
    end


%% Plot marker positions using raw data
figure
    for posture = 1:2
        for target = 1:8
            traj = trialData(min(find(visualLabel==2 & proprioLabel==posture & directionLabel==target))).kinematics_updated.marker.position;
            if posture == 1
                plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
            else
                plot(traj(:,1),traj(:,2),'--','Color',tcmap(target,:));
            end
            hold on
        end
    end

%% Plot cursor positions using raw data
figure
    for posture = 1:2
        for target = 1:8
            traj = trialData(min(find(visualLabel==2 & proprioLabel==posture & directionLabel==target))).kinematics_updated.cursor.position;
            if posture == 1
                plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
            else
                plot(traj(:,1),traj(:,2),'--','Color',tcmap(target,:));
            end
            hold on
        end
    end

%% Plot Positions using trajStruct (processed data)
figure
    for posture = 1:2
        for target = 1:8
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj;
            if posture == 1
                plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
            else
                plot(traj(:,1),traj(:,2),'--','Color',tcmap(target,:));
            end
            hold on
        end
    end
    
