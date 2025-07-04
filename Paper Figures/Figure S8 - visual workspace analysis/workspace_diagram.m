clear; clc;

%% Load and unpack raw data file
    dataset = 'R20200221_all_visual_data';
    rawData = load('D:\Animals\Rocky\2020\02\20200221\Rocky_20200221_preprocessedData_20220512_165624.mat');
    trialData = [rawData.trialData];
    visualLabel = [trialData.visualLabel];
    proprioLabel = [trialData.proprioLabel];
    directionLabel = [trialData.directionLabel];

%% Get trajStruct for makrer trajectories
    [Data,zScoreParams] = loadData(dataset);
    
    %Get trajStruct
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    trajFields = {'zSmoothFR','marker'};
    trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);                     
    

%% Setup colormap (based on number of postures)
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
    [pcmap,tcmap,rainbow] = getColorMaps(numPostures);        
    
%% Get cursor and target radius
    cursorRadius = rawData.taskInfo.cursorRadius;
    targetRadius = rawData.taskInfo.targetRadius;
    
%% Get target locations for all task settings
    targetLocationStruct = struct('posture',[],'visualID',[],'targetLocations',[]);
    structInd = 1;
    
    for posture = 1:2
        for visualID = 1:2
            targetLocations = zeros(8,2);
            for target = 1:8
                tempData = trialData(min(find(visualLabel==visualID & proprioLabel==posture & directionLabel==target))).z_extraInfo.endTargetPosition(1:2).*1000;  
                targetLocations(target,:) = tempData;
            end
            targetLocationStruct(structInd).posture = posture;
            targetLocationStruct(structInd).visualID = visualID;
            targetLocationStruct(structInd).targetLocations = targetLocations;
            structInd = structInd + 1;
        end
    end
    
%% Plot average marker trajectories for all task settings

    figure
    hold on
    for posture = 1:2
        for visualID = 1:2
            for target = 1:8
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target & [trajStruct.visual]==visualID).avgMarker.traj;
                if visualID == 1
                    plot(traj(:,1),traj(:,2),'Color',pcmap(posture,:));
                else
                    plot(traj(:,1),traj(:,2),'--','Color',pcmap(posture,:));
                end
            end
        end
    end
    
    
%% Plot all target locations (from raw data)
    %figure;
    hold on;
    for posture = 1:2
        for visualID = 1:2
            cursorTransform = trialData(min(find(visualLabel==visualID & proprioLabel==posture))).z_extraInfo.cursorTransform(1:2,4)*1000;
            
            targetLocations = targetLocationStruct([targetLocationStruct.visualID]==visualID & [targetLocationStruct.posture]==posture).targetLocations;
            for target = 1:8
                circle(targetLocations(target,1)-cursorTransform(1), targetLocations(target,2)-cursorTransform(2), targetRadius, 'k');
            end
        end
    end
    axis('equal')

    
%% Plot larger targets to make figure
    figure;
    figureTargetRadius = 10;
    hold on;
    for posture = 1:2
        for visualID = 1:2
            cursorTransform = trialData(min(find(visualLabel==visualID & proprioLabel==posture))).z_extraInfo.cursorTransform(1:2,4)*1000;
            
            targetLocations = targetLocationStruct([targetLocationStruct.visualID]==visualID & [targetLocationStruct.posture]==posture).targetLocations;
            for target = 1:8
                circle(targetLocations(target,1)-cursorTransform(1), targetLocations(target,2)-cursorTransform(2), figureTargetRadius, 'k');
            end
        end
    end
    axis('equal')
