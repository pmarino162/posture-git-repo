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

%     for posture = 1:2
%         posture
%         cursorTransform = trialData(min(find(visualLabel==2 & proprioLabel==posture))).z_extraInfo.cursorTransform(1:2,4)*1000
%         for target = 1:8
%             tempData = trialData(min(find(visualLabel==2 & proprioLabel==posture & directionLabel==target))).z_extraInfo.endTargetPosition(1:2).*1000;
%             tempData - cursorTransform'
%         end
%     end

% %% Plot all marker trajectories
% % 
%     % Get overall workspace center & dist bt workspaces
%     workspaceCenter1 = trialData(min(find(visualLabel==2 & proprioLabel==1))).z_extraInfo.cursorTransform(1:2,4)*1000;
%     workspaceCenter2 = trialData(min(find(visualLabel==2 & proprioLabel==2))).z_extraInfo.cursorTransform(1:2,4)*1000;
%     xRange = abs(workspaceCenter1(1)-workspaceCenter2(1));
%     yRange = abs(workspaceCenter1(2)-workspaceCenter2(2));
%     xMin = min(workspaceCenter1(1),workspaceCenter2(1));
%     yMin = min(workspaceCenter1(2),workspaceCenter2(2));
%     workspaceCenter = [xMin+xRange/2;yMin+yRange/2];
%     workspaceCenter1 = workspaceCenter1 - workspaceCenter;
%     workspaceCenter2 = workspaceCenter2 - workspaceCenter;
%     workspaceDist = norm(workspaceCenter1-workspaceCenter2);
%     
%     % Get target distance
%     startTarget = trialData(1).z_extraInfo.startTargetPosition;
%     endTarget = trialData(1).z_extraInfo.endTargetPosition;
%     targetDist = norm(endTarget-startTarget)*1000;
%     
% 
% 
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
%     %% Print target location in marker coordinates
%     for posture = 1:2
%         posture
%         cursorTransform = trialData(min(find(visualLabel==2 & proprioLabel==posture))).z_extraInfo.cursorTransform(1:2,4)*1000
%         for target = 1:8
%             tempData = trialData(min(find(visualLabel==2 & proprioLabel==posture & directionLabel==target))).z_extraInfo.endTargetPosition(1:2).*1000;
%             tempData - cursorTransform'
%         end
%     end
% 
% 
% %% Plot marker positions using raw data
% figure
%     for posture = 1:2
%         for target = 1:8
%             traj = trialData(min(find(visualLabel==2 & proprioLabel==posture & directionLabel==target))).kinematics_updated.marker.position;
%             if posture == 1
%                 plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
%             else
%                 plot(traj(:,1),traj(:,2),'--','Color',tcmap(target,:));
%             end
%             hold on
%         end
%     end
% 
% %% Plot cursor positions using raw data
% figure
%     for posture = 1:2
%         for target = 1:8
%             traj = trialData(min(find(visualLabel==2 & proprioLabel==posture & directionLabel==target))).kinematics_updated.cursor.position;
%             if posture == 1
%                 plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
%             else
%                 plot(traj(:,1),traj(:,2),'--','Color',tcmap(target,:));
%             end
%             hold on
%         end
%     end

