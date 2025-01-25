clear; clc; clf; close all;

%% Setup saveFig    saveFig = true;
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22';
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = customRainbow;
    tAngleCmap = interp1([1;100],[1 1 1; 1 0 0],[1:1:100]');
    
%% Load Data, subselect, get traj struct
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;


    task = 'reaching';
    [Data] = loadEarlData20210706_20211210; 
    
    %Keep only trials with delay length >= 500ms
    kinData = [Data.kinData];
    delayLength = [kinData.delayLength];
    rmTrials = delayLength < 500;
    Data(rmTrials) = [];
    
    %Reach time histogram (Exclude +/-1.5*std)
    clearvars rmTrials
    kinData = [Data.kinData];
    reachTime = [kinData.reachTime];
    meanReachTime = mean(reachTime);
    stdReachTime = std(reachTime);
    figure; histogram(reachTime);
    hold on
    ax = gca;
    line([meanReachTime meanReachTime],[ax.YLim],'LineWidth',2,'Color','b')
    hold on
    acceptDist = 1.5*stdReachTime;
    line([meanReachTime+acceptDist meanReachTime+acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
    line([meanReachTime-acceptDist meanReachTime-acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
    numTrials = size(Data,2);
    rmTrials = reachTime < meanReachTime-acceptDist | reachTime > meanReachTime+acceptDist; 
    numTrialsRemaining = sum(numTrials - sum(rmTrials));
    text(0.95*ax.XLim(2),0.95*ax.YLim(2),[num2str(numTrialsRemaining),'/',num2str(numTrials),...
        ' (',num2str(round(100*numTrialsRemaining/numTrials)),'%) ',' remaining'],'HorizontalAlignment','right');
    xlabel('Reach Time (ms)')
    ylabel('Count')
    %Exclude trials based on reach time
    Data(rmTrials) = [];
    
    %Examine remaining condition counts
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialCounts = countConditionTrials(Data,condFields);
    numCondTrials = [trialCounts.numTrials];
    figure
    histogram(numCondTrials)
    xlabel('Number of Trials')
    ylabel('Number of Conditions')

    %Remove conditions based on condition counts
    rmCond = [trialCounts.numTrials] < 6;
    rmTrialNums = [trialCounts(rmCond).trialNum];
    trialNum = [Data.trialNum];
    Data(ismember(trialNum,rmTrialNums)) = [];
    close all
    
    trialInclStates(1).trialName = {'GridReaching'};
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%     trajStruct = trajStruct(ismember([trajStruct.posture],[1:4]));

   trialInclStates(1).inclStates = {{'state','Delay','first',250},{'state','Target Acquire','first',0}};
   delayTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%  %%trajStruct = trajStruct(ismember([trajStruct.posture],[1:4]));

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

%% Get number of postures; max and min postural differences
    switch task
        %Reaching and Planning
        case {'reaching','planning'}
            postureList = unique([trajStruct.posture]);
        postureData = [Data.conditionData];
        postureData  = [postureData.postureID];
        workspaceCenter = [Data.targetData];
        workspaceCenter = vertcat(workspaceCenter.workspaceCenter);
        postureDiffStruct = struct('posture1',[],'posture2',[],'postureDiff',[]);
        structInd = 1;
        for posture1Ind = 1:length(postureList)
            posture1 = postureList(posture1Ind);
            p1WC = workspaceCenter(find(postureData==posture1,1),:);
            for posture2Ind = 1:length(postureList)
                posture2 = postureList(posture2Ind);
                p2WC = workspaceCenter(find(postureData==posture2,1),:);
                postureDiff = vecnorm(p1WC-p2WC);
                postureDiffStruct(structInd).posture1 = posture1;
                postureDiffStruct(structInd).posture2 = posture2;
                postureDiffStruct(structInd).postureDiff = postureDiff;
                structInd = structInd + 1;
            end
        end
        postureDiff = [postureDiffStruct.postureDiff];
        postureDiff = postureDiff(postureDiff > 0);
%         minPostureDiff = min(postureDiff);
        minPostureDiff = 0;
        maxPostureDiff = max(postureDiff);
    end
    

%% Do PCA on condition averages for visualization; get CIs
    avgSmoothFR = [trajStruct.avgSmoothFR];
    allAvgs = vertcat(avgSmoothFR.traj);
    [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);       
    % Add Projections to trajStruct 
    for i = 1:size(trajStruct,2)
       trajStruct(i).PCA = (trajStruct(i).avgSmoothFR.traj-mu)*coeff;
       PCA = trajStruct(i).PCA;
       numTrials = size(trajStruct(i).allSmoothFR,2);
       allPCA = NaN(minNumTimestamps,size(PCA,2),numTrials);
       for j = 1:numTrials
          trajStruct(i).allPCA(j).traj = (trajStruct(i).allSmoothFR(j).traj-mu)*coeff;
          trajStruct(i).allPCA(j).timestamps = trajStruct(i).allSmoothFR(j).timestamps;
          allPCA(:,:,j) = trajStruct(i).allPCA(j).traj(1:minNumTimestamps,:);
       end
       %Get confidence intervals
       trajStruct(i).PCACI = 1.96*std(allPCA,0,3)/sqrt(numTrials);
    end    
    
%% Do the heavy lifting
    % Create distance distributions
    resultStruct = struct('cond1P',[],'cond1T',[],'cond2P',[],'cond2T',[],'dist',[],'alpha',[],'fullAlpha',[],'uDist',[],'uAlpha',[],'uFullAlpha',[],'targetAngle',[],'postureDiff',[]);
    numCond = size(trajStruct,2);
    numDraws = 1; %Number of random draws from each condition
    numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw
    numPts = minNumTimestamps; %Number of points from each trajectory
    
    structInd = 1;
    for cond1 = 1:numCond
        cond1
        for cond2 = cond1:numCond
            cond1P = trajStruct(cond1).posture;
            cond1T = trajStruct(cond1).target;
            cond2P = trajStruct(cond2).posture;
            cond2T = trajStruct(cond2).target;
            condStr = num2str([cond1P,cond1T,cond2P,cond2T]);
            resultStruct(structInd).cond1P = cond1P;
            resultStruct(structInd).cond1T = cond1T;
            resultStruct(structInd).cond2P = cond2P;
            resultStruct(structInd).cond2T = cond2T;
            if strcmpi(task,'BCI') | strcmpi(task,'iso') | strcmpi(task,'multijoint BCI')
                resultStruct(structInd).postureDiff = abs(cond2P-cond1P);
            elseif strcmpi(task,'reaching') | strcmpi(task,'planning') 
                resultStruct(structInd).postureDiff = postureDiffStruct([postureDiffStruct.posture1]==cond1P & [postureDiffStruct.posture2]==cond2P).postureDiff;
            end
            targetAngle = max([cond1T,cond2T])-min([cond1T,cond2T]);
            if targetAngle > 4
                targetAngle = 8-targetAngle;
            end
            resultStruct(structInd).targetAngle = round(100*targetAngle./4);
            for i = 1:numDraws
                %Create traj1
                numTraj1 = size(trajStruct(cond1).allSmoothFR,2);
                sampInd1 = randsample(numTraj1,numSample);
                traj1Struct = trajStruct(cond1).allSmoothFR(sampInd1);
                traj1 = getAvgTraj20211210(traj1Struct,binWidth);
                %Create traj2
                if cond1 == cond2
                    numTrajRemaining = numTraj1 - numSample;
                    sampInd2 = randsample(numTrajRemaining,numSample);
                    remainingInd = setdiff(1:numTraj1,sampInd1);
                    sampInd2 = remainingInd(sampInd2);
                else
                    numTraj2 = size(trajStruct(cond2).allSmoothFR,2);
                    sampInd2 = randsample(numTraj2,numSample);
                end
                traj2Struct = trajStruct(cond2).allSmoothFR(sampInd2);
                traj2 = getAvgTraj20211210(traj2Struct,binWidth);
                %Get traj1 and traj2 mean prep state
                trialNumTraj1 = [trajStruct(cond1).allSmoothFR(sampInd1).trialNum];
                trialNumTraj2 = [trajStruct(cond2).allSmoothFR(sampInd2).trialNum];
                prepStatesTraj1 = NaN(numSample,size(traj1,2));
                prepStatesTraj2 = NaN(numSample,size(traj2,2));
                for sample = 1:numSample
                   delayCond1TrajInd = find([delayTrajStruct(cond1).allSmoothFR.trialNum]==trialNumTraj1(sample)); 
                   prepStatesTraj1(sample,:) = mean(delayTrajStruct(cond1).allSmoothFR(delayCond1TrajInd).traj,1); 
                   delayCond2TrajInd = find([delayTrajStruct(cond2).allSmoothFR.trialNum]==trialNumTraj2(sample)); 
                   prepStatesTraj2(sample,:) = mean(delayTrajStruct(cond2).allSmoothFR(delayCond2TrajInd).traj,1); 
                end
                alpha = mean(prepStatesTraj1,1) - mean(prepStatesTraj2,1);
                %Shift traj2 by alpha
                unshiftTraj2 = traj2;
                traj2 = traj2 + alpha;
                resultStruct(structInd).fullAlpha(i,:) = alpha;
                alpha = vecnorm(alpha);
                resultStruct(structInd).alpha(i) = alpha;
                %Get residual distance
                dist = getMeanDist(traj1,traj2,numPts);
                resultStruct(structInd).dist(i) = dist;
                %Visualize trajectories
                visualizeTraj(traj1,traj2,mu,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist,alpha)
            end
            resultStruct(structInd).uDist = mean([resultStruct(structInd).dist]);
            resultStruct(structInd).uAlpha = mean([resultStruct(structInd).alpha]);
            resultStruct(structInd).uFullAlpha = mean([resultStruct(structInd).fullAlpha]);
            structInd = structInd + 1;
        end
    end
    
%% Collect lists of condition variables
    cond1P = [resultStruct.cond1P]; cond1T = [resultStruct.cond1T]; cond2P = [resultStruct.cond2P]; cond2T = [resultStruct.cond2T];

%% Make target/posture comparison plot
    %Within condition comparisons
    withinConditionResult = resultStruct([cond1P == cond2P] & [cond1T == cond2T]);
    withinConditionDist = [withinConditionResult.uDist];
    withinConditionAlpha = [withinConditionResult.uAlpha];
    withinConditionPostureDiff = [withinConditionResult.postureDiff];
    withinConditionTargetAngle = [withinConditionResult.targetAngle];
    withinConditionTargetAngle = withinConditionTargetAngle + 1;
    %Across posture comparisons
    acrossPostureResult = resultStruct([cond1P ~= cond2P] & [cond1T == cond2T]);
    acrossPostureDist = [withinConditionDist,[acrossPostureResult.uDist]];
    acrossPostureAlpha = [withinConditionAlpha,[acrossPostureResult.uAlpha]];
    acrossPostureDiff = [withinConditionPostureDiff,[acrossPostureResult.postureDiff]];
    rPosture = corrcoef([acrossPostureAlpha',acrossPostureDist']);
    rPosture = rPosture(1,2);
    XPosture = [ones(length(acrossPostureAlpha),1),acrossPostureAlpha'];
    YPosture = acrossPostureDist';
    BPosture = inv(XPosture'*XPosture)*XPosture'*YPosture;

    %Across target comparisons
    acrossTargetResult = resultStruct([cond1P == cond2P] & [cond1T ~= cond2T]);
    acrossTargetDist = [withinConditionDist,[acrossTargetResult.uDist]];
    acrossTargetAlpha = [withinConditionAlpha,[acrossTargetResult.uAlpha]];
    acrossTargetAngle = [withinConditionTargetAngle,[acrossTargetResult.targetAngle]];
    rTarget = corrcoef([acrossTargetAlpha',acrossTargetDist']);
    rTarget = rTarget(1,2);
    XTarget = [ones(length(acrossTargetAlpha),1),acrossTargetAlpha'];
    YTarget = acrossTargetDist';
    BTarget = inv(XTarget'*XTarget)*XTarget'*YTarget;

    f = figure; f.Position = [40 40 1000 750];
    hold on
    for i = 1:length(acrossTargetAlpha)
        plot(acrossTargetAlpha(i),acrossTargetDist(i),'o','MarkerSize',10,...
            'MarkerEdgeColor','r','MarkerFaceColor',tAngleCmap(acrossTargetAngle(i),:))
    end
    for i = 1:length(acrossPostureAlpha)
        postureDiff = acrossPostureDiff(i);
        pDiffColor = interp1([minPostureDiff;maxPostureDiff],[1 1 1; 0 0 1],postureDiff);
        plot(acrossPostureAlpha(i),acrossPostureDist(i),'o','MarkerSize',10,...
            'MarkerEdgeColor','b','MarkerFaceColor',pDiffColor);
    end
    ax=gca;
    xLim = ax.XLim;
    temp = [ones(2,1),xLim'];
    plot(xLim,temp*BPosture,'b','LineWidth',1)
    plot(xLim,temp*BTarget,'r','LineWidth',1)
    xlabel('Shift (Hz)')
    ylabel('Trajectory Warping (Hz)')
    yLim = ax.YLim;
    yRange = yLim(2)-yLim(1); 

  
    if saveFig
       saveas(gcf,fullfile(saveDir,[task,' Figures'],'ShiftDistComparisonWithText.jpg')); 
    end
    
    

%% Define local functions
    %Get optimal alpha
    function alpha = getOptimalAlpha(traj1,traj2,numPts)
        numDims = size(traj1,2);
        alpha = NaN(1,numDims);
        for dim = 1:numDims
            alpha(dim) = (1/numPts)*sum(traj1(1:numPts,dim)-traj2(1:numPts,dim));
        end
    end
    
    %Get mean dist
    function dist = getMeanDist(traj1,traj2,numPts)
        dist = 0;
        for i = 1:numPts
            dist = vecnorm(traj1(i,:)-traj2(i,:)) + dist;
        end
        dist = dist/numPts;
    end
    
    %Visualize traj
    function [] = visualizeTraj(traj1,traj2,mu,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist,alpha)
        traj1Proj = (traj1-mu)*coeff;
        traj2Proj = (traj2-mu)*coeff;
        unshiftTraj2Proj = (unshiftTraj2-mu)*coeff;
        traj1AvgProj = trajStruct(cond1).PCA;
        traj2AvgProj = trajStruct(cond2).PCA;
        
        figure
        hold on
        %Traj 1
        plot3(traj1Proj(:,1),traj1Proj(:,2),traj1Proj(:,3),'Color','r');
        for pt = 1:numPts
            ptInd = mod(pt,8);
            if ptInd == 0
                ptInd = 8;
            end
            plot3(traj1Proj(pt,1),traj1Proj(pt,2),traj1Proj(pt,3),'.','MarkerSize',10,'Color',tcmap(ptInd,:));
        end
        %Unshifted Traj 2
        plot3(unshiftTraj2Proj(:,1),unshiftTraj2Proj(:,2),unshiftTraj2Proj(:,3),'Color','b');
        %Shifted Traj 2
        plot3(traj2Proj(:,1),traj2Proj(:,2),traj2Proj(:,3),'Color','g');
        for pt = 1:numPts
            ptInd = mod(pt,8);
            if ptInd == 0
                ptInd = 8;
            end
            plot3(traj2Proj(pt,1),traj2Proj(pt,2),traj2Proj(pt,3),'.','MarkerSize',10,'Color',tcmap(ptInd,:));
        end
        %Avg Traj 1
        plot3(traj1AvgProj(:,1),traj1AvgProj(:,2),traj1AvgProj(:,3),'Color','r','LineWidth',2);
        plot3(traj1AvgProj(1,1),traj1AvgProj(1,2),traj1AvgProj(1,3),'.','MarkerSize',10,'Color','r');
        %Avg Traj 2
        plot3(traj2AvgProj(:,1),traj2AvgProj(:,2),traj2AvgProj(:,3),'Color','b','LineWidth',2);
        plot3(traj2AvgProj(1,1),traj2AvgProj(1,2),traj2AvgProj(1,3),'.','MarkerSize',10,'Color','b');

        %Get condition info 
        cond1P = trajStruct(cond1).posture;
        cond1T = trajStruct(cond1).target;
        cond2P = trajStruct(cond2).posture;
        cond2T = trajStruct(cond2).target;
        
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
        title(['P',num2str(cond1P),'T',num2str(cond1T),' vs ','P',num2str(cond2P),'T',num2str(cond2T)...
            newline,'Dist = ',num2str(dist),'Hz',newline,'Alpha = ',num2str(alpha)])

        close 
    end