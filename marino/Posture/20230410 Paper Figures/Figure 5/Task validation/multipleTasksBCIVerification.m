clear; clc; clf; close all;

%% Colormap
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    
%% Load data; get trajStruct
    %Load data
    dataset = 'E20200314';
    [Data,zScoreParams] = loadData(dataset);
    %Keep only BCI
    conditionData = [Data.conditionData]; taskID = [conditionData.taskID];
    Data = Data([taskID==1]);
    %Get trajStruct
    binWidth = 25; kernelStdDev = 25;
    trialInclStates = struct('trialName','','inclStates',[]);
    trajFields = {'zSmoothFR','bciCursorTraj'};
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    
    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
    trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',0},{'state','Step 2','first',0}};

    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

%% For each posture, get all target locations
    targetStruct = struct('posture',[],'centerLoc',[],'centerSize',[],'targetLoc',[],'workspaceCenter',[]);
    postureList = unique([trajStruct.posture]);
    conditionData = [Data.conditionData];
    postureData = [conditionData.postureID];
    structInd = 1;
    for posture = postureList
        tempData = Data(postureData==posture);
        targetData = [tempData.targetData];
        targetStruct(structInd).posture = posture;
        targetStruct(structInd).centerLoc = unique(vertcat(targetData.centerLoc),'rows');
        targetStruct(structInd).centerSize = unique(vertcat(targetData.centerSize),'rows');
        
        targetIDs = [targetData.targetID];
        targetList = unique(targetIDs);
        targetLoc = struct('targetID',[],'targetLoc',[],'targetSize',[]);
        targetLocInd = 1;
        for targetID = targetList
            targetLoc(targetLocInd).targetID = targetID;
            tempTargetData = targetData(targetIDs==targetID);
            targetLoc(targetLocInd).targetLoc = unique(vertcat(tempTargetData.target1Loc),'rows');
            targetLoc(targetLocInd).targetSize = unique(vertcat(tempTargetData.target1Size),'rows');
            targetLocInd = targetLocInd + 1;
        end
        
        targetStruct(structInd).targetLoc = targetLoc;
        targetStruct(structInd).workspaceCenter = unique(vertcat(targetData.workspaceCenter),'rows');
        structInd = structInd + 1;
    end

%% Plot kinematic traces
    %plot traces
    targetList = unique([trajStruct.target]);    
    for posture = postureList
        figure; hold on;
        title(['Posture ',num2str(posture)]);
        xlabel('x (mm)')
        ylabel('y (mm)');
        
        for target = targetList
            %Plot individual trials
            allTraj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allBciCursorTraj;
            numTrials = size(allTraj,2);
%             for trial = 1:numTrials
%                 traj = allTraj(trial).traj;
%                 plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
%             end
            %Plot average
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgBciCursorTraj.traj;
            plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
        end
        
        %add target info
        for target = targetList
            targetLocStruct = targetStruct([targetStruct.posture]==posture).targetLoc;
            targetLoc = targetLocStruct([targetLocStruct.targetID]==target).targetLoc;
            targetSize = targetLocStruct([targetLocStruct.targetID]==target).targetSize;
            circle(targetLoc(:,1),targetLoc(:,2),targetSize,'k');

        end
        %add center info
        centerLoc = targetStruct([targetStruct.posture]==posture).centerLoc;
        centerSize = targetStruct([targetStruct.posture]==posture).centerSize;     
        circle(centerLoc(:,1),centerLoc(:,2),centerSize,'k');
        
    end
    

%% 