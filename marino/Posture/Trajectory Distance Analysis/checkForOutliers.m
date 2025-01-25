clear; clc; clf; close all;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = customRainbow;
    pAngleCmap = interp1([25;100],[1 1 1; 0 0 1],[1:1:100]');
    tAngleCmap = interp1([25;100],[1 1 1; 1 0 0],[1:1:100]');
    
%% Load Data, subselect, get traj struct
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;

    % BCI
    [Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
    [Data] = subselectForTrajDist(Data,'BCI');
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-250},{'state','Step 2','first',0}};
    trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%     %Across posture BCI
%     [Data] = loadEarlData20210901_20211210;
%     [Data] = subselectForTrajDist(Data,'Multi Joint BCI');
    
    
%       %Iso force
%     [Data] = loadEarlData20200116_20211210();
%     [Data] = subselectForTrajDist(Data,'iso');
%     trialInclStates(1).trialName = {'IsometricForce_1D'};
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-100},{'state','Target Hold','first',0}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%     
% %     Reaching
%     [Data] = loadEarlData20210706_20211210; 
%     [Data] = subselectForTrajDist(Data,'reaching');
%     trialInclStates(1).trialName = {'GridReaching'};
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-150},{'state','Target Hold','first',0}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%     trajStruct = trajStruct(ismember([trajStruct.posture],[1:4]));
    %     Planning
%     [Data] = loadEarlData20210706_20211210; 
%     [Data] = subselectForTrajDist(Data,'planning');
%      trialInclStates(1).trialName = {'GridReaching'};
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'state','Delay','first',0},{'state','Delay','first',250}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%  trajStruct = trajStruct(ismember([trajStruct.posture],[1:4]));

%% Get timestamps dist and min number timestamps
    numTimestamps = [];
    numCondTraj = [];
    for i = 1:size(trajStruct,2)
       numTraj = size(trajStruct(i).allSmoothFR,2);
       numCondTraj = [numCondTraj,numTraj];
       for j = 1:numTraj
          numTimestamps = [numTimestamps,size(trajStruct(i).allSmoothFR(j).timestamps,2)]; 
       end
    end
    figure
    histogram(numTimestamps)
    xlabel('Number of 25ms bins')
    ylabel('Number of trials')
    
    figure
    histogram(numCondTraj)
    xlabel('Number of trials')
    ylabel('Number of conditions')
    
% Get minimum number of trials and timestamps
    [minNumTimestamps,i] = min(numTimestamps);
    [minNumCondTraj,i] = min(numCondTraj);


%% Do PCA on condition averages for visualization 
    avgSmoothFR = [trajStruct.avgSmoothFR];
    allAvgs = vertcat(avgSmoothFR.traj);
    [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);       
    % Add Projections to trajStruct 
    for i = 1:size(trajStruct,2)
       trajStruct(i).PCA = trajStruct(i).avgSmoothFR.traj*coeff;
       numTrials = size(trajStruct(i).allSmoothFR,2);
       for j = 1:numTrials
          trajStruct(i).allPCA(j).traj = trajStruct(i).allSmoothFR(j).traj*coeff;
          trajStruct(i).allPCA(j).timestamps = trajStruct(i).allSmoothFR(j).timestamps;
       end
    end    
    
%% Plot PCA to look for outliers
posture = 2;
target = 7;

allPCA = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allPCA;
numTrials = size(allPCA,2);

figure
for trial = 1:numTrials
    traj = allPCA(trial).traj;
    time = allPCA(trial).timestamps;
    for i = 1:10
        subplot(2,5,i)
        plot(time,traj(:,i))
        hold on
    end
    
end

%% Compare two condition averages
cond1P = 1;
cond1T = 7;
cond2P = 2;
cond2T = 7;

cond1Traj = trajStruct([trajStruct.posture]==cond1P & [trajStruct.target]==cond1T).PCA;
cond1Time = trajStruct([trajStruct.posture]==cond1P & [trajStruct.target]==cond1T).avgSmoothFR.timestamps;
cond2Traj = trajStruct([trajStruct.posture]==cond2P & [trajStruct.target]==cond2T).PCA;
cond2Time = trajStruct([trajStruct.posture]==cond2P & [trajStruct.target]==cond2T).avgSmoothFR.timestamps;

figure
for i = 1:10
    subplot(2,5,i)
    plot(cond1Time,cond1Traj(:,i));
    hold on
    plot(cond2Time,cond2Traj(:,i));
    
end

%Get axis ranges
    for plotInd = 1:10
       ax = subplot(2,5,plotInd);
       range(plotInd) = ax.YLim(1,2)-ax.YLim(1,1);
       center(plotInd) = (ax.YLim(1,2)+ax.YLim(1,1))/2;
    end
    maxRange = max(range);
    for plotInd = 1:10
       ax = subplot(2,5,plotInd);
       ax.YLim = [center(plotInd)-maxRange/2 center(plotInd)+maxRange/2];
       xlabel('time (ms)')
        ylabel(['PC ',num2str(plotInd)])
    end
