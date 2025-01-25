clear; clc; clf; close all;

%% Load data; get trajStruct
    %Load data
    dataset = 'E20200314';
    [Data,zScoreParams] = loadData(dataset);
    %Keep only iso force
    conditionData = [Data.conditionData]; taskID = [conditionData.taskID];
    Data = Data([taskID==3]);
    %Get trajStruct
    binWidth = 25; kernelStdDev = 25;
    trialInclStates = struct('trialName','','inclStates',[]);
    trajFields = {'force','forceCursor'};
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'IsometricForce_1D'};
    trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Success with Reward','first',0}};
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);


%% For each posture, get all target locations
    targetStruct = struct('posture',[],'centerLoc',[],'centerSize',[],'targetLoc',[],...
        'workspaceCenter',[]);
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
            targetLoc(targetLocInd).targetLoc = unique(vertcat(tempTargetData.targetLoc),'rows');
            targetLoc(targetLocInd).targetSize = unique(vertcat(tempTargetData.targetSize),'rows');
            targetLocInd = targetLocInd + 1;
        end
        
        targetStruct(structInd).targetLoc = targetLoc;
        targetStruct(structInd).workspaceCenter = unique(vertcat(targetData.workspaceCenter),'rows');
        structInd = structInd + 1;
    end

%% Plot force v time with target locations
    %plot traces
    targetList = unique([trajStruct.target]);    
    for posture = postureList
        figure; hold on;
        title(['Posture ',num2str(posture)]);
        xlabel('time (ms)')
        ylabel('F_Y (N)');
        
        for target = targetList
            %Plot individual trials
            allTraj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allForceCursor;
            numTrials = size(allTraj,2);
            for trial = 1:numTrials
                traj = allTraj(trial).traj(:,2);
                time = allTraj(trial).timestamps;
                if target==3
                    plot(time,traj,'r-')
                elseif target == 7
                    plot(time,traj,'b-')
                end
            end
            %Plot average
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgForceCursor.traj(:,2);
            time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgForceCursor.timestamps;
            if target==3
                plot(time,traj,'r-','LineWidth',2)
            elseif target == 7
                plot(time,traj,'b-','LineWidth',2)
            end
        end
        
        %add target info
        ax = gca;
        plotXLims = ax.XLim; 
        for target = targetList
            targetLocStruct = targetStruct([targetStruct.posture]==posture).targetLoc;
            targetLoc = targetLocStruct([targetLocStruct.targetID]==target).targetLoc(:,2);
            targetSize = targetLocStruct([targetLocStruct.targetID]==target).targetSize(:,2);
            plot(plotXLims,(targetLoc+targetSize/2)*ones(1,2),'k');
            plot(plotXLims,(targetLoc-targetSize/2)*ones(1,2),'k');
        end
        %add center info
        centerLoc = targetStruct([targetStruct.posture]==posture).centerLoc(:,2);
        centerSize = targetStruct([targetStruct.posture]==posture).centerSize(:,2);     
        plot(plotXLims,(centerLoc+centerSize/2)*ones(1,2),'k');
        plot(plotXLims,(centerLoc-centerSize/2)*ones(1,2),'k');
        
    end
    
