clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\Consistent Posture Signal';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;    
    
%% Load Data, subselect
    dataset = 'E20200314';
    [Data,zScoreParams] = loadData(dataset);
    
    %Get trajStruct
    binWidth = 25; kernelStdDev = 25;
    trialInclStates = struct('trialName','','inclStates',[]);
    trajFields = {'zSmoothFR','markerVel','marker'};
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    
    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
        trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',25},{'state','Step 2','first',0}};
        %trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',25},{'state','Step 2','first',150}};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-300},{'kin','moveOnsetTime','first',-50}};
        %trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-300},{'state','Hold','first',0}};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {{'state','Target','first',75},{'state','Target','first',200}};
        %trialInclStates(3).inclStates = {{'state','Target','first',25},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

    taskIDs = struct('ID',[],'task','');
    taskIDs(1).ID = 1; taskIDs(1).task = 'BC';
    taskIDs(2).ID = 2; taskIDs(2).task = 'HC';
    taskIDs(3).ID = 3; taskIDs(3).task = 'Iso';
    
%% Eliminate conditions with too few trials; get minNumCondTraj
    numCondTraj = [];
    for i = 1:size(trajStruct,2)
       numTraj = size(trajStruct(i).allSmoothFR,2);
       numCondTraj = [numCondTraj,numTraj];
    end

    figure
    histogram(numCondTraj)
    xlabel('Number of trials')
    ylabel('Number of conditions')
    
    %Get rid of conditions with less than 22 trials
    minNumCondTraj = 22;
    trajStruct = trajStruct(numCondTraj >= minNumCondTraj);
    
%% Do PCA on condition averages for visualization
    avgSmoothFR = [trajStruct.avgSmoothFR];
    allAvgs = vertcat(avgSmoothFR.traj);
    [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);       
    % Add Projections to trajStruct 
    for i = 1:size(trajStruct,2)
       trajStruct(i).PCA = (trajStruct(i).avgSmoothFR.traj)*coeff;
       PCA = trajStruct(i).PCA;
       numTrials = size(trajStruct(i).allSmoothFR,2);
       for j = 1:numTrials
          trajStruct(i).allPCA(j).traj = (trajStruct(i).allSmoothFR(j).traj-mu)*coeff;
          trajStruct(i).allPCA(j).timestamps = trajStruct(i).allSmoothFR(j).timestamps;
       end
    end   
    
%% For each task, get compositionality indices and lower bound
    resultStruct = struct('cond1Task',[],'cond2Task',[],'cond1P',[],'cond1T',[],'cond2P',[],'cond2T',[],'dist',[],'alpha',[],'fullAlpha',[],'uDist',[],'uAlpha',[],'uFullAlpha',[],'targetAngle',[],'postureDiff',[]);
    numCond = size(trajStruct,2);
    numDraws = 50; %Number of random draws from each condition
    numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw
    structInd = 1;
    
    for cond1 = 1:numCond
        cond1
        for cond2 = cond1:numCond
            cond1Task = trajStruct(cond1).task;
            cond1P = trajStruct(cond1).posture;
            cond1T = trajStruct(cond1).target;
            cond2Task = trajStruct(cond2).task;
            cond2P = trajStruct(cond2).posture;
            cond2T = trajStruct(cond2).target;
            condStr = num2str([cond1Task,cond2Task,cond1P,cond1T,cond2P,cond2T]);
            resultStruct(structInd).cond1Task = cond1Task;
            resultStruct(structInd).cond1P = cond1P;
            resultStruct(structInd).cond1T = cond1T;
            resultStruct(structInd).cond2Task = cond2Task;
            resultStruct(structInd).cond2P = cond2P;
            resultStruct(structInd).cond2T = cond2T;
            resultStruct(structInd).postureDiff = abs(cond2P-cond1P);
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
                %Get numPts = minNumTimestamps bt 2 trajectories
                numPts = min(size(traj1,1),size(traj2,1));
                %Get alpha
                %alpha = getOptimalAlpha(traj1,traj2,numPts);
                alpha = traj1(1,:)-traj2(1,:);
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
                %visualizeTraj(traj1,traj2,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist,alpha)
            end
            resultStruct(structInd).uDist = mean([resultStruct(structInd).dist]);
            resultStruct(structInd).uAlpha = mean([resultStruct(structInd).alpha]);
            resultStruct(structInd).uFullAlpha = mean([resultStruct(structInd).fullAlpha]);
            structInd = structInd + 1;
        end
    end
    
%% Collect lists of condition variables
    cond1Task = [resultStruct.cond1Task]; cond1P = [resultStruct.cond1P]; cond1T = [resultStruct.cond1T]; cond2Task = [resultStruct.cond2Task]; cond2P = [resultStruct.cond2P]; cond2T = [resultStruct.cond2T];    

%% Plot results
    %Within condition comparisons
    withinConditionResult = resultStruct([cond1Task == cond2Task] & [cond1P == cond2P] & [cond1T == cond2T]);
    withinConditionDist = [withinConditionResult.uDist];
    withinConditionAlpha = [withinConditionResult.uAlpha];
    withinConditionPostureDiff = [withinConditionResult.postureDiff];
    withinConditionTargetAngle = [withinConditionResult.targetAngle];
    withinConditionTargetAngle = withinConditionTargetAngle + 1;

    %Across posture comparisons (within Task)
    acrossPostureResult = resultStruct([cond1Task == cond2Task] & [cond1P ~= cond2P] & [cond1T == cond2T]);
    acrossPostureDist = [[acrossPostureResult.uDist]];
    acrossPostureAlpha = [[acrossPostureResult.uAlpha]];
    acrossPostureDiff = [[acrossPostureResult.postureDiff]];
    
    %Across target comparisons (within Task)
    acrossTargetResult = resultStruct([cond1Task == cond2Task] & [cond1P == cond2P] & [cond1T ~= cond2T]);
    acrossTargetDist = [[acrossTargetResult.uDist]];
    acrossTargetAlpha = [[acrossTargetResult.uAlpha]];
    acrossTargetDiff = [[acrossTargetResult.postureDiff]];    
    
    bciResult = resultStruct([cond1Task == 1] & [cond2Task == 1] & [cond1T == cond2T]);
    bciAlpha = [bciResult.uAlpha];
    bciDist = [bciResult.uDist];
    bciX = [ones(length(bciAlpha),1),bciAlpha'];
    bciY = bciDist';
    bciB = inv(bciX'*bciX)*bciX'*bciY;
    
    isoResult = resultStruct([cond1Task == 3] & [cond2Task == 3] & [cond1T == cond2T]);
    isoAlpha = [isoResult.uAlpha];
    isoDist = [isoResult.uDist];
    isoX = [ones(length(isoAlpha),1),isoAlpha'];
    isoY = isoDist';
    isoB = inv(isoX'*isoX)*isoX'*isoY;
    
    reachResult = resultStruct([cond1Task == 2] & [cond2Task == 2] & [cond1T == cond2T]);
    reachAlpha = [reachResult.uAlpha];
    reachDist = [reachResult.uDist];
    reachX = [ones(length(reachAlpha),1),reachAlpha'];
    reachY = reachDist';
    reachB = inv(reachX'*reachX)*reachX'*reachY;
    
    
    
    f = figure; f.Position = [40 40 1000 750];
    hold on
    taskCmap = [1 0 0; 0 1 0; 0 0 1];
    for i = 1:length(withinConditionResult)
           task = withinConditionResult(i).cond1Task;
           plot(withinConditionAlpha(i),withinConditionDist(i),'o','MarkerSize',10,...
            'MarkerEdgeColor','g','MarkerFaceColor',taskCmap(task,:));
    end
    
    for i = 1:length(acrossPostureResult)
           task = acrossPostureResult(i).cond1Task;
           plot(acrossPostureAlpha(i),acrossPostureDist(i),'o','MarkerSize',10,...
            'MarkerEdgeColor','b','MarkerFaceColor',taskCmap(task,:));
    end
    
%     for i = 1:length(acrossTargetResult)
%            task = acrossTargetResult(i).cond1Task;
%            plot(acrossTargetAlpha(i),acrossTargetDist(i),'s','MarkerSize',10,...
%             'MarkerEdgeColor','b','MarkerFaceColor',taskCmap(task,:));
%     end

    ax=gca;
    xLim = ax.XLim;
    temp = [ones(2,1),xLim'];
    plot(xLim,temp*bciB,'r','LineWidth',1)
     plot(xLim,temp*isoB,'b','LineWidth',1)
     plot(xLim,temp*reachB,'g','LineWidth',1)
    yLim = ax.YLim;
    yRange = yLim(2)-yLim(1); 

xlabel('Postural Shift')
ylabel('Trajectory Warping')

%% Make     
% %% Get stats
%     bciCompInd = resultStruct([resultStruct.task]==1).compInd; 
%     reachCompInd = resultStruct([resultStruct.task]==2).compInd; 
%     isoCompInd = resultStruct([resultStruct.task]==3).compInd; 
%     [h,p] = ttest2(bciCompInd,reachCompInd,'Vartype','unequal')
%     [h,p] = ttest2(bciCompInd,isoCompInd,'Vartype','unequal')
%     
% %% Plot
%     figure
%     hold on
%     for task = 1:3
%        compInd = resultStruct([resultStruct.task]==task).compInd;
%        plot(task*ones(size(compInd)),compInd,'.','MarkerSize',10);
%        plot(task,mean(compInd),'.','MarkerSize',30) 
%     end
%     xlim([0,4])
%     xticks([1:3])
%     xticklabels({'BCI','Reaching','IsoForce'})
    
%% Local functions 
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
    function [] = visualizeTraj(traj1,traj2,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist,alpha)
        traj1Proj = traj1*coeff;
        traj2Proj = traj2*coeff;
        unshiftTraj2Proj = unshiftTraj2*coeff;
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
        cond1Task = trajStruct(cond1).task;
        cond1P = trajStruct(cond1).posture;
        cond1T = trajStruct(cond1).target;
        cond2Task = trajStruct(cond2).task;
        cond2P = trajStruct(cond2).posture;
        cond2T = trajStruct(cond2).target;
        
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
        title(['Task',num2str(cond1Task),'P',num2str(cond1P),'T',num2str(cond1T),' vs ','Task',num2str(cond2Task),'P',num2str(cond2P),'T',num2str(cond2T)...
            newline,'Dist = ',num2str(dist),'Hz',newline,'Alpha = ',num2str(alpha)])

        close 
    end
%% Archived code

%  for task = unique([trajStruct.task])
%         %Get task traj struct
%         taskTrajStruct = trajStruct([trajStruct.task]==task);
%         %Do PCA on task trajectories for visualization later
%         %[coeff,taskTrajStruct] = pcaForTaskVis(taskTrajStruct);
%         %Main Loop
%         compInd = [];
%         taskTargetList = unique([taskTrajStruct.target]);
%         for target = taskTargetList
%            taskTargetPostureList = unique([taskTrajStruct([taskTrajStruct.target]==target).posture]);
%            for posture1 = taskTargetPostureList
%               posture1Ind = find(taskTargetPostureList==posture1);
%               for posture2 = taskTargetPostureList(posture1Ind:end)
%                   compIndDraws = [];
%                   for draw = 1:numDraws
%                       %Create traj1
%                       traj1Struct = taskTrajStruct([taskTrajStruct.target]==target & [taskTrajStruct.posture]==posture1);
%                       numTraj1 = size(traj1Struct.allSmoothFR,2);
%                       sampInd1 = randsample(numTraj1,numSample);
%                       traj1 = getAvgTraj20211210(traj1Struct.allSmoothFR(sampInd1),binWidth);
%                       %Create traj2
%                       traj2Struct = taskTrajStruct([taskTrajStruct.target]==target & [taskTrajStruct.posture]==posture2);
%                       numTraj2 = size(traj2Struct.allSmoothFR,2);
%                       sampInd2 = randsample(numTraj2,numSample);
%                       traj2 = getAvgTraj20211210(traj2Struct.allSmoothFR(sampInd2),binWidth);
%                       %Get minNumTimestamps bt the 2 trajectories for comparison
%                       minNumTimestamps = min(size(traj1,1),size(traj2,1));
%                       %Get posture shift
%                       postureShift = traj1(1,:)-traj2(1,:);
%                       postureDist = vecnorm(postureShift);
%                       %Shift traj2 by posture shift
%                       traj2 = traj2 + postureShift;
%                       %Get residual distance
%                       dist = getMeanDist(traj1,traj2,minNumTimestamps);
%                       curCompInd = dist/abs(posture2-posture1);
%                       %Visualize sampled trajectories
%                       %visualizeTraj(traj1,traj2,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist,alpha)
%                       %Save results
%                       compIndDraws = [compIndDraws,curCompInd];
%                   end
%                   %Get mean across draws; store
%                   if posture1 == posture2
%                       lowerBound = [lowerBound,mean(compIndDraws)];
%                   else
%                       compInd = [compInd,mean(compIndDraws)];
%                   end
%               end
%            end
%         end
%         resultStruct(structInd).task = task;
%         resultStruct(structInd).compInd = compInd;
%         structInd = structInd + 1;
%     end


% %PCA for visualization 
% function [coeff,taskTrajStruct] = pcaForTaskVis(taskTrajStruct)
%     avgSmoothFR = [taskTrajStruct.avgSmoothFR];
%     allAvgs = vertcat(avgSmoothFR.traj);
%     [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);       
%     % Add Projections to taskTrajStruct 
%     for i = 1:size(taskTrajStruct,2)
%        taskTrajStruct(i).PCA = (taskTrajStruct(i).avgSmoothFR.traj)*coeff;
%        PCA = taskTrajStruct(i).PCA;
%        numTrials = size(taskTrajStruct(i).allSmoothFR,2);
%        allPCA = NaN(minNumTimestamps,size(PCA,2),numTrials);
%        for j = 1:numTrials
%           taskTrajStruct(i).allPCA(j).traj = (taskTrajStruct(i).allSmoothFR(j).traj-mu)*coeff;
%           taskTrajStruct(i).allPCA(j).timestamps = taskTrajStruct(i).allSmoothFR(j).timestamps;
%           allPCA(:,:,j) = taskTrajStruct(i).allPCA(j).traj(1:minNumTimestamps,:);
%        end
%     end   
% end