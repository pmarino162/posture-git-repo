clear; clc; clf; close all;

%% Load, Preprocess, and Label Data    
    [Data] = loadEarlData20210828;
    Data = Data([Data.trialStatus]==1);
    
%% Create trajStruct
    condFields = {{'posture','conditionData','postureID'},{'task','conditionData','taskID'},{'target','targetData','targetID'}};
    trajFields = {'marker'};
    
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {0};
        trialInclStates(1).addTimeToEnd = {0};
%     trialInclStates(2).trialName = {'GridReaching'};
    trialInclStates(2).trialName = {'DelayedCenterOut20210828'};
        trialInclStates(2).inclStates = {'No Cheating','Target Acquire'};
        trialInclStates(2).inclOccurrence = {'last','last'};
        trialInclStates(2).addTimeToBeginning = {0,0};
        trialInclStates(2).addTimeToEnd = {0,0};    
        
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
            
%% Collect Data
    tempData = trajStruct([trajStruct.posture]==2 & [trajStruct.task]==1);
    N00_BCHandPose = [];
    for i=1:size(tempData,2)
        N00_BCHandPose = vertcat(N00_BCHandPose,tempData(i).avgMarker.traj);
    end
    
    tempData = trajStruct([trajStruct.posture]==1 & [trajStruct.task]==1);
    I30_BCHandPose = [];
    for i=1:size(tempData,2)
        I30_BCHandPose = vertcat(I30_BCHandPose,tempData(i).avgMarker.traj);
    end
    
    tempData = trajStruct([trajStruct.posture]==3 & [trajStruct.task]==1);
    E30_BCHandPose = [];
    for i=1:size(tempData,2)
        E30_BCHandPose = vertcat(E30_BCHandPose,tempData(i).avgMarker.traj);
    end
    
    tempData = trajStruct([trajStruct.posture]==2 & [trajStruct.task]==2);
    N00_HCHandPose = [];
    for i=1:size(tempData,2)
        N00_HCHandPose = vertcat(N00_HCHandPose,tempData(i).avgMarker.traj);
    end
    
    

%% Load Colormaps
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    cmap = cmap([1,3,5],:);
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;

%% Plot Data
    figure
    hcCenterLoc = [-40,-420];
    
    HCTrajStruct = trajStruct([trajStruct.task]==2);
    targetList = [HCTrajStruct.target];
    for target = targetList
        traj = HCTrajStruct([HCTrajStruct.target]==target).avgMarker.traj;
        plot(traj(:,1),traj(:,2),'.','MarkerSize',10,'Color',tcmap(target,:));
        hold on;
    end
    
    
    tempData = N00_BCHandPose;
    for i = 1:size(tempData,1)
        plot(tempData(i,1)-hcCenterLoc(1),tempData(i,2)-hcCenterLoc(2),'.','MarkerSize',20,'Color',cmap(2,:));
        hold on;
    end
    
    tempData = I30_BCHandPose;
    for i = 1:size(tempData,1)
        plot(tempData(i,1)-hcCenterLoc(1),tempData(i,2)-hcCenterLoc(2),'.','MarkerSize',20,'Color',cmap(1,:));
        hold on;
    end

    tempData = E30_BCHandPose;
    for i = 1:size(tempData,1)
        plot(tempData(i,1)-hcCenterLoc(1),tempData(i,2)-hcCenterLoc(2),'.','MarkerSize',20,'Color',cmap(3,:));
        hold on;
    end
    
axis equal
%     axis square
    xlabel('x (mm)')
    ylabel('y (mm)')

%% Print Means