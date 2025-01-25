clear; clc; clf; close all;

%% Setup saveFig    saveFig = true;
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220126 - traj difference vs shift';
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = customRainbow;
    tAngleCmap = interp1([25;100],[1 1 1; 1 0 0],[1:1:100]');
    
%% Load Data, subselect, get traj struct
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;

%     % BCI
%     task = 'BCI';
%     [Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
%     [Data] = subselectForTrajDist(Data,task);
%     trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
%     condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
% %     trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-150},{'state','Step 2','first',0}};
%     trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%     %Across posture BCI
%     [Data] = loadEarlData20210901_20211210;
%     [Data] = subselectForTrajDist(Data,'Multi Joint BCI');
    
    
%       %Iso force
%     task = 'iso';
%     [Data] = loadEarlData20200116_20211210();
%     [Data] = subselectForTrajDist(Data,task);
%     trialInclStates(1).trialName = {'IsometricForce_1D'};
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%     
%     Reaching
        task = 'reaching';
    [Data] = loadEarlData20210706_20211210; 
    [Data] = subselectForTrajDist(Data,task);
    trialInclStates(1).trialName = {'GridReaching'};
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%     trajStruct = trajStruct(ismember([trajStruct.posture],[1:4]));
%     trajStruct = trajStruct(ismember([trajStruct.posture],[1,2,3,5,6,7]));


% % %         Planning
% %     task = 'planning';
% %     [Data] = loadEarlData20210706_20211210; 
% %     [Data] = subselectForTrajDist(Data,task);
% %      trialInclStates(1).trialName = {'GridReaching'};
% %     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
% %     trialInclStates(1).inclStates = {{'state','Delay','first',50},{'state','Delay','first',250}};
% %     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
% %  %%trajStruct = trajStruct(ismember([trajStruct.posture],[1:4]));

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
        %BCI and Iso
        case {'BCI','iso'}
        maxPostureDiff = 4;
        minPostureDiff = 1;
        
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
        minPostureDiff = min(postureDiff);
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
    numDraws = 100; %Number of random draws from each condition
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
            if strcmpi(task,'BCI') | strcmpi(task,'iso') 
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
                %Get alpha
                alpha = getOptimalAlpha(traj1,traj2,numPts);
%                 alpha = traj1(1,:)-traj2(1,:);
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
%                 visualizeTraj(traj1,traj2,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist,alpha)
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
    rWC = corrcoef([withinConditionAlpha',withinConditionDist']);
    rWC = rWC(1,2);
    XWC = [ones(length(withinConditionAlpha),1),withinConditionAlpha'];
    YWC = withinConditionDist';
    BWC = inv(XWC'*XWC)*XWC'*YWC;

    %Across posture comparisons
    acrossPostureResult = resultStruct([cond1P ~= cond2P] & [cond1T == cond2T]);
    acrossPostureDist = [acrossPostureResult.uDist];
    acrossPostureAlpha = [acrossPostureResult.uAlpha];
    acrossPostureDiff = [acrossPostureResult.postureDiff];
    rPosture = corrcoef([acrossPostureAlpha',acrossPostureDist']);
    rPosture = rPosture(1,2);
    XPosture = [ones(length(acrossPostureAlpha),1),acrossPostureAlpha'];
    YPosture = acrossPostureDist';
    BPosture = inv(XPosture'*XPosture)*XPosture'*YPosture;

    %Across target comparisons
    acrossTargetResult = resultStruct([cond1P == cond2P] & [cond1T ~= cond2T]);
    acrossTargetDist = [acrossTargetResult.uDist];
    acrossTargetAlpha = [acrossTargetResult.uAlpha];
    acrossTargetAngle = [acrossTargetResult.targetAngle];
    rTarget = corrcoef([acrossTargetAlpha',acrossTargetDist']);
    rTarget = rTarget(1,2);
    XTarget = [ones(length(acrossTargetAlpha),1),acrossTargetAlpha'];
    YTarget = acrossTargetDist';
    BTarget = inv(XTarget'*XTarget)*XTarget'*YTarget;

    f = figure; f.Position = [40 40 1400 900];
        subplot(5,5,1)
        plot(acrossTargetAlpha,acrossTargetDist,'.','Color','r')
        hold on
        ax=gca;
        xLim = ax.XLim;
        temp = [ones(2,1),xLim'];
        plot(xLim,temp*BTarget,'k','LineWidth',2)
        xlabel('Shift (Hz)')
        ylabel('Trajectory Difference (Hz)')
        title(['Change Target',newline,'r = ',num2str(round(rTarget,2)),...
            ' W = ',num2str(round(BTarget(2),2))])
    subplot(5,5,11)
        plot(acrossPostureAlpha,acrossPostureDist,'.','Color','b')
        hold on
        ax=gca;
        xLim = ax.XLim;
        temp = [ones(2,1),xLim'];
        plot(xLim,temp*BPosture,'k','LineWidth',2)
        xlabel('Shift (Hz)')
        ylabel('Trajectory Difference (Hz)')
        title(['Change Posture',newline,'r = ',num2str(round(rPosture,2)),...
            ' W = ',num2str(round(BPosture(2),2))])
    subplot(5,5,21)
        plot(withinConditionAlpha,withinConditionDist,'.','Color','g')
        hold on
        ax=gca;
        xLim = ax.XLim;
        temp = [ones(2,1),xLim'];
        plot(xLim,temp*BWC,'k','LineWidth',2)
        xlabel('Shift (Hz)')
        ylabel('Trajectory Difference (Hz)')
        title(['Within Condition',newline,'r = ',num2str(round(rWC,2)),...
            ' W = ',num2str(round(BWC(2),2))])
    subplot(5,5,[3 4 5 8 9 10 13 14 15 18 19 20]) 
        hold on
        for i = 1:length(withinConditionAlpha)
           	plot(withinConditionAlpha(i),withinConditionDist(i),'o','MarkerSize',10,...
                'MarkerEdgeColor','g','MarkerFaceColor',[1 1 1]);
        end
        for i = 1:length(acrossTargetAlpha)
            plot(acrossTargetAlpha(i),acrossTargetDist(i),'o','MarkerSize',10,...
                'MarkerEdgeColor','r','MarkerFaceColor',tAngleCmap(acrossTargetAngle(i),:))
            condStr = [acrossTargetResult(i).cond1P,acrossTargetResult(i).cond1T,acrossTargetResult(i).cond2P,acrossTargetResult(i).cond2T];
            text(acrossTargetAlpha(i),acrossTargetDist(i),num2str(condStr));
        end
        for i = 1:length(acrossPostureAlpha)
            postureDiff = acrossPostureDiff(i);
            pDiffColor = interp1([minPostureDiff;maxPostureDiff],[1 1 1; 0 0 1],postureDiff);
            plot(acrossPostureAlpha(i),acrossPostureDist(i),'o','MarkerSize',10,...
                'MarkerEdgeColor','b','MarkerFaceColor',pDiffColor);
            condStr = [acrossPostureResult(i).cond1P,acrossPostureResult(i).cond1T,acrossPostureResult(i).cond2P,acrossPostureResult(i).cond2T];
            text(acrossPostureAlpha(i),acrossPostureDist(i),num2str(condStr));
        end
        ax=gca;
        xLim = ax.XLim;
        temp = [ones(2,1),xLim'];
        plot(xLim,temp*BWC,'g','LineWidth',1)
        plot(xLim,temp*BPosture,'b','LineWidth',1)
        plot(xLim,temp*BTarget,'r','LineWidth',1)
        xlabel('Shift (Hz)')
        ylabel('Trajectory Difference (Hz)')
        yLim = ax.YLim;
        yRange = yLim(2)-yLim(1); 
   subplot(5,5,[2 7 12 17])
        binWidth = 2;
        histogram(withinConditionDist,[yLim(1):binWidth:yLim(2)],'Normalization','probability','FaceColor','g')
        hold on
        histogram(acrossPostureDist,[yLim(1):binWidth:yLim(2)],'Normalization','probability','FaceColor','b')
        histogram(acrossTargetDist,[yLim(1):binWidth:yLim(2)],'Normalization','probability','FaceColor','r')
        view([270 90])
        ax=gca;
        ax.XLim = yLim;
        xticklabels([])
        yticklabels([])     
  subplot(5,5,[23 24 25])
        binWidth = 2;
        histogram(withinConditionAlpha,[xLim(1):binWidth:xLim(2)],'Normalization','probability','FaceColor','g')
        hold on
        histogram(acrossPostureAlpha,[xLim(1):binWidth:xLim(2)],'Normalization','probability','FaceColor','b')
        histogram(acrossTargetAlpha,[xLim(1):binWidth:xLim(2)],'Normalization','probability','FaceColor','r')
        xticklabels([])
        yticklabels([])
  
    if saveFig
       saveas(gcf,fullfile(saveDir,task,'ShiftDistComparisonWithText.fig')); 
    end
    
    
     f = figure; f.Position = [40 40 1400 900];
        subplot(5,5,1)
        plot(acrossTargetAlpha,acrossTargetDist,'.','Color','r')
        hold on
        ax=gca;
        xLim = ax.XLim;
        temp = [ones(2,1),xLim'];
        plot(xLim,temp*BTarget,'k','LineWidth',2)
        xlabel('Shift (Hz)')
        ylabel('Trajectory Difference (Hz)')
        title(['Change Target',newline,'r = ',num2str(round(rTarget,2)),...
            ' W = ',num2str(round(BTarget(2),2))])
    subplot(5,5,11)
        plot(acrossPostureAlpha,acrossPostureDist,'.','Color','b')
        hold on
        ax=gca;
        xLim = ax.XLim;
        temp = [ones(2,1),xLim'];
        plot(xLim,temp*BPosture,'k','LineWidth',2)
        xlabel('Shift (Hz)')
        ylabel('Trajectory Difference (Hz)')
        title(['Change Posture',newline,'r = ',num2str(round(rPosture,2)),...
            ' W = ',num2str(round(BPosture(2),2))])
    subplot(5,5,21)
        plot(withinConditionAlpha,withinConditionDist,'.','Color','g')
        hold on
        ax=gca;
        xLim = ax.XLim;
        temp = [ones(2,1),xLim'];
        plot(xLim,temp*BWC,'k','LineWidth',2)
        xlabel('Shift (Hz)')
        ylabel('Trajectory Difference (Hz)')
        title(['Within Condition',newline,'r = ',num2str(round(rWC,2)),...
            ' W = ',num2str(round(BWC(2),2))])
    subplot(5,5,[3 4 5 8 9 10 13 14 15 18 19 20]) 
        hold on
        for i = 1:length(withinConditionAlpha)
           	plot(withinConditionAlpha(i),withinConditionDist(i),'o','MarkerSize',10,...
                'MarkerEdgeColor','g','MarkerFaceColor',[1 1 1]);
        end
        for i = 1:length(acrossTargetAlpha)
            plot(acrossTargetAlpha(i),acrossTargetDist(i),'o','MarkerSize',10,...
                'MarkerEdgeColor','r','MarkerFaceColor',tAngleCmap(acrossTargetAngle(i),:))
            condStr = [acrossTargetResult(i).cond1P,acrossTargetResult(i).cond1T,acrossTargetResult(i).cond2P,acrossTargetResult(i).cond2T];
        end
        for i = 1:length(acrossPostureAlpha)
            postureDiff = acrossPostureDiff(i);
            pDiffColor = interp1([minPostureDiff;maxPostureDiff],[1 1 1; 0 0 1],postureDiff);
            plot(acrossPostureAlpha(i),acrossPostureDist(i),'o','MarkerSize',10,...
                'MarkerEdgeColor','b','MarkerFaceColor',pDiffColor);
            condStr = [acrossPostureResult(i).cond1P,acrossPostureResult(i).cond1T,acrossPostureResult(i).cond2P,acrossPostureResult(i).cond2T];
        end
        ax=gca;
        xLim = ax.XLim;
        temp = [ones(2,1),xLim'];
        plot(xLim,temp*BWC,'g','LineWidth',1)
        plot(xLim,temp*BPosture,'b','LineWidth',1)
        plot(xLim,temp*BTarget,'r','LineWidth',1)
        xlabel('Shift (Hz)')
        ylabel('Trajectory Difference (Hz)')
        yLim = ax.YLim;
        yRange = yLim(2)-yLim(1); 
   subplot(5,5,[2 7 12 17])
        binWidth = 2;
        histogram(withinConditionDist,[yLim(1):binWidth:yLim(2)],'Normalization','probability','FaceColor','g')
        hold on
        histogram(acrossPostureDist,[yLim(1):binWidth:yLim(2)],'Normalization','probability','FaceColor','b')
        histogram(acrossTargetDist,[yLim(1):binWidth:yLim(2)],'Normalization','probability','FaceColor','r')
        view([270 90])
        ax=gca;
        ax.XLim = yLim;
        xticklabels([])
        yticklabels([])     
  subplot(5,5,[23 24 25])
        binWidth = 2;
        histogram(withinConditionAlpha,[xLim(1):binWidth:xLim(2)],'Normalization','probability','FaceColor','g')
        hold on
        histogram(acrossPostureAlpha,[xLim(1):binWidth:xLim(2)],'Normalization','probability','FaceColor','b')
        histogram(acrossTargetAlpha,[xLim(1):binWidth:xLim(2)],'Normalization','probability','FaceColor','r')
        xticklabels([])
        yticklabels([])
  
    if saveFig
       saveas(gcf,fullfile(saveDir,task,'ShiftDistComparison.jpg')); 
    end
    close all
%% Visualize some averages
    cond1P = 2;
    cond1T = 7;
    cond2P = 5;
    cond2T = 7;

    cond1Traj = trajStruct([trajStruct.posture]==cond1P & [trajStruct.target]==cond1T).PCA(1:minNumTimestamps,:);
    cond1Time = trajStruct([trajStruct.posture]==cond1P & [trajStruct.target]==cond1T).avgSmoothFR.timestamps(1:minNumTimestamps);
    cond1CI =  trajStruct([trajStruct.posture]==cond1P & [trajStruct.target]==cond1T).PCACI;
    cond2Traj = trajStruct([trajStruct.posture]==cond2P & [trajStruct.target]==cond2T).PCA(1:minNumTimestamps,:);
    cond2Time = trajStruct([trajStruct.posture]==cond2P & [trajStruct.target]==cond2T).avgSmoothFR.timestamps(1:minNumTimestamps);
    cond2CI = trajStruct([trajStruct.posture]==cond2P & [trajStruct.target]==cond2T).PCACI;
    dist = resultStruct([resultStruct.cond1P]==cond1P & [resultStruct.cond1T]==cond1T & [resultStruct.cond2P]==cond2P & [resultStruct.cond2T]==cond2T).uDist;
    alpha = resultStruct([resultStruct.cond1P]==cond1P & [resultStruct.cond1T]==cond1T & [resultStruct.cond2P]==cond2P & [resultStruct.cond2T]==cond2T).uAlpha;
    fullAlpha = resultStruct([resultStruct.cond1P]==cond1P & [resultStruct.cond1T]==cond1T & [resultStruct.cond2P]==cond2P & [resultStruct.cond2T]==cond2T).uFullAlpha;
    cond2ShiftTraj = trajStruct([trajStruct.posture]==cond2P & [trajStruct.target]==cond2T).avgSmoothFR.traj(1:minNumTimestamps,:);
    cond2ShiftTraj = cond2ShiftTraj + fullAlpha;
    cond2ShiftTraj = (cond2ShiftTraj-mu)*coeff;

    %Unshifted
    f = figure; f.Position = [100 100 1500 500];
    for i = 1:10
        subplot(2,5,i)
        shadedErrorBar(cond1Time,cond1Traj(:,i),cond1CI(:,i),'lineprops',{'LineWidth',2})
        hold on
        shadedErrorBar(cond2Time,cond2Traj(:,i),cond2CI(:,i),'lineprops',{'LineWidth',2})
    end
    sgtitle(['P',num2str(cond1P),'T',num2str(cond1T),' vs ','P',num2str(cond2P),'T',num2str(cond2T),' Unshifted',...
        newline,'Dist = ',num2str(round(dist)),' Shift = ',num2str(round(alpha))])
    legend({['P',num2str(cond1P),'T',num2str(cond1T)],['P',num2str(cond2P),'T',num2str(cond2T)]})
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
        
    if saveFig
       saveas(gcf,fullfile(saveDir,task,['P',num2str(cond1P),'T',num2str(cond1T),'vs','P',num2str(cond2P),'T',num2str(cond2T),'TrajComparisonUnshifted.jpg'])); 
    end
%     close all
    
    
    
    %Shifted
    f = figure; f.Position = [100 100 1500 500];
    for i = 1:10
        subplot(2,5,i)
        shadedErrorBar(cond1Time,cond1Traj(:,i),cond1CI(:,i),'lineprops',{'LineWidth',2})
        hold on
        shadedErrorBar(cond2Time,cond2ShiftTraj(:,i),cond2CI(:,i),'lineprops',{'LineWidth',2})
    end
    sgtitle(['P',num2str(cond1P),'T',num2str(cond1T),' vs ','P',num2str(cond2P),'T',num2str(cond2T),' Shifted',...
        newline,'Dist = ',num2str(round(dist)),' Shift = ',num2str(round(alpha))])
    legend({['P',num2str(cond1P),'T',num2str(cond1T)],['P',num2str(cond2P),'T',num2str(cond2T)]})
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
    if saveFig
       saveas(gcf,fullfile(saveDir,task,['P',num2str(cond1P),'T',num2str(cond1T),'vs','P',num2str(cond2P),'T',num2str(cond2T),'TrajComparisonShifted.jpg'])); 
    end
%     close all
    
%% Look for outliers within a condition        
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