clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\ReachOrthQuant';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    
%% Load data; get trajstruct
    dataset = 'E20210706';
    [Data,zScoreParams] = loadData(dataset);
    centerHoldData = Data;
    targetHoldData = Data;
        
    %Get trajStruct
    binWidth = 50;
    kernelStdDev = 50;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    %All
    trialInclStates(1).trialName = {'GridReaching'};
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',500}};
    %Reach
    reachTrialInclStates(1).trialName = {'GridReaching'};
    reachCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    reachTrialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
    %Center Hold
    centerHoldTrialInclStates(1).trialName = {'GridReaching'};
    centerHoldCondFields = {{'posture','conditionData','postureID'}};
    centerHoldTrialInclStates(1).inclStates = {{'state','Center Hold','first',100},{'state','Delay','first',50}};
    %Target Hold
    targetHoldTrialInclStates(1).trialName = {'GridReaching'};
    targetHoldCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    targetHoldTrialInclStates(1).inclStates = {{'state','Target Hold','first',100},{'state','Success with Reward','first',50}};
    %Target Hold Plot 
    targetHoldPlotTrialInclStates(1).trialName = {'GridReaching'};
    targetHoldPlotCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    targetHoldPlotTrialInclStates(1).inclStates = {{'state','Target Hold','first',-200},{'state','Success with Reward','first',0}};
    %Get TrajStruct
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    reachTrajStruct = getTrajStruct20220419(Data,reachCondFields,trajFields,reachTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    centerHoldTrajStruct = getTrajStruct20220419(centerHoldData,centerHoldCondFields,trajFields,centerHoldTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    targetHoldTrajStruct = getTrajStruct20220419(targetHoldData,targetHoldCondFields,trajFields,targetHoldTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    targetHoldPlotTrajStruct = getTrajStruct20220419(targetHoldData,targetHoldPlotCondFields,trajFields,targetHoldPlotTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    
%% Keep only Postures 1-4
    trajStruct = trajStruct(ismember([trajStruct.posture],1:4));
    reachTrajStruct = reachTrajStruct(ismember([reachTrajStruct.posture],1:4));
    centerHoldTrajStruct = centerHoldTrajStruct(ismember([centerHoldTrajStruct.posture],1:4));
    targetHoldTrajStruct = targetHoldTrajStruct(ismember([targetHoldTrajStruct.posture],1:4));
    targetHoldPlotTrajStruct = targetHoldPlotTrajStruct(ismember([targetHoldPlotTrajStruct.posture],1:4));
    
%% Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);
    
%% Find Center Hold Posture Axis 
    trajSplitLim = round(750/binWidth);
    % For each timecourse in trajstruct, get single point
    for i = 1:size(centerHoldTrajStruct,2)
        %Individual trials
        for j = 1:size(centerHoldTrajStruct(i).allSmoothFR,2)
           traj = centerHoldTrajStruct(i).allSmoothFR(j).traj;
           trajLength = size(traj,1);
           if trajLength >= trajSplitLim
               splitPt = round(trajLength/2);
               centerHoldTrajStruct(i).allSmoothFR(j).trialAvg1 = mean(centerHoldTrajStruct(i).allSmoothFR(j).traj(1:splitPt,:)); 
               centerHoldTrajStruct(i).allSmoothFR(j).trialAvg2 = mean(centerHoldTrajStruct(i).allSmoothFR(j).traj(splitPt:end,:)); 
           else
               centerHoldTrajStruct(i).allSmoothFR(j).trialAvg1 = mean(centerHoldTrajStruct(i).allSmoothFR(j).traj); 
           end
        end
    end
    for i = 1:size(targetHoldTrajStruct,2)
        %Individual trials
        for j = 1:size(targetHoldTrajStruct(i).allSmoothFR,2)
           traj = targetHoldTrajStruct(i).allSmoothFR(j).traj;
           trajLength = size(traj,1);
           if trajLength >= trajSplitLim
               splitPt = round(trajLength/2);
               targetHoldTrajStruct(i).allSmoothFR(j).trialAvg1 = mean(targetHoldTrajStruct(i).allSmoothFR(j).traj(1:splitPt,:)); 
               targetHoldTrajStruct(i).allSmoothFR(j).trialAvg2 = mean(targetHoldTrajStruct(i).allSmoothFR(j).traj(splitPt:end,:)); 
           else
               targetHoldTrajStruct(i).allSmoothFR(j).trialAvg1 = mean(targetHoldTrajStruct(i).allSmoothFR(j).traj); 
           end
        end
    end
    
    % For each posture, get all pts of each type
    allObsStruct = struct('posture',[],'allObsCH',[],'numObsCH',[],'allObsTH',[],'numObsTH',[]);
    structInd = 1;
    for posture = postureList
        %Center hold obs
        allObsCH = [];
        tempTrajStruct = centerHoldTrajStruct([centerHoldTrajStruct.posture]==posture);
        for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               trialAvg1 = tempTrajStruct(i).allSmoothFR(j).trialAvg1;
               trialAvg2 = tempTrajStruct(i).allSmoothFR(j).trialAvg2;
               allObsCH = vertcat(allObsCH,vertcat(trialAvg1,trialAvg2));
           end
        end
        %Get Target Hold Data (Leftward reaches)
        allObsTH = [];
        targetHoldPosture = posture + 1;
        targetHoldTarget = 5;
        tempTrajStruct = targetHoldTrajStruct([targetHoldTrajStruct.posture]==targetHoldPosture & [targetHoldTrajStruct.target]==targetHoldTarget);
        for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               trialAvg1 = tempTrajStruct(i).allSmoothFR(j).trialAvg1;
               trialAvg2 = tempTrajStruct(i).allSmoothFR(j).trialAvg2;
               allObsTH = vertcat(allObsTH,vertcat(trialAvg1,trialAvg2));
           end
        end
       %Get Target Hold Data (Rightward reaches)
       targetHoldPosture = posture -1;
       targetHoldTarget = 1;
       tempTrajStruct = targetHoldTrajStruct([targetHoldTrajStruct.posture]==targetHoldPosture & [targetHoldTrajStruct.target]==targetHoldTarget);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               trialAvg1 = tempTrajStruct(i).allSmoothFR(j).trialAvg1;
               trialAvg2 = tempTrajStruct(i).allSmoothFR(j).trialAvg2;
               allObsTH = vertcat(allObsTH,vertcat(trialAvg1,trialAvg2));
           end
       end
       allObsStruct(structInd).posture = posture;
       allObsStruct(structInd).allObsCH = allObsCH;
       allObsStruct(structInd).numObsCH = size(allObsCH,1);
       allObsStruct(structInd).allObsTH = allObsTH;
       allObsStruct(structInd).numObsTH = size(allObsTH,1);
       structInd = structInd + 1;
    end
    
    %Create obsStruct using balanced points
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;
    
    for posture = postureList
        allObs = [];
        allObsCH = allObsStruct([allObsStruct.posture]==posture).allObsCH;
        numObsCH = allObsStruct([allObsStruct.posture]==posture).numObsCH;
        allObsTH = allObsStruct([allObsStruct.posture]==posture).allObsTH;
        numObsTH = allObsStruct([allObsStruct.posture]==posture).numObsTH;
        minNumSamples = min(size(allObsCH,1),size(allObsTH,1));
        allObs = vertcat(allObs,datasample(allObsCH,minNumSamples,'Replace',false));
        allObs = vertcat(allObs,datasample(allObsTH,minNumSamples,'Replace',false));
        obsStruct(structInd).label = posture;
        obsStruct(structInd).allObs = allObs;
        obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
        structInd = structInd + 1;
    end
    
    %Do PCA on observations
    allObs = vertcat(obsStruct.allObs);
    [coeff,score,latent,tsquared,explained,mu] = pca(allObs);
    
    %Replace obsStruct observations w PC projections
    numPCs = 20;
    for i = 1:numel(obsStruct)
       obsStruct(i).allObs = obsStruct(i).allObs*coeff(:,1:numPCs);       
    end
    
    %Do LDA on obsStruct
    [postureLDA] = doLDA(obsStruct);
    
    %Combine PCA and LDA to get posture Axis 
    postureLDA = coeff(:,1:numPCs)*postureLDA;
    
    % Do LDA by posture on all trial averages
%     obsStruct = struct('label',[],'numObs',[],'allObs',[]);
%     structInd = 1;   
%     for posture = postureList
%        obsStruct(structInd).label = posture;
%        allObs = []; allObsCH = []; allObsTH = [];
%        %Get Center Hold Data
%        tempTrajStruct = centerHoldTrajStruct([centerHoldTrajStruct.posture]==posture);
%        for i = 1:size(tempTrajStruct,2)
%            for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
%                trialAvg = tempTrajStruct(i).allSmoothFR(j).trialAvg;
%                allObsCH = vertcat(allObsCH,trialAvg);
%            end
%        end
%        %Get Target Hold Data (Leftward reaches)
%        targetHoldPosture = posture + 1;
%        targetHoldTarget = 5;
%        tempTrajStruct = targetHoldTrajStruct([targetHoldTrajStruct.posture]==targetHoldPosture & [targetHoldTrajStruct.target]==targetHoldTarget);
%        for i = 1:size(tempTrajStruct,2)
%            for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
%                trialAvg = tempTrajStruct(i).allSmoothFR(j).trialAvg;
%                allObsTH = vertcat(allObsTH,trialAvg);
%            end
%        end
%        %Get Target Hold Data (Rightward reaches)
%        targetHoldPosture = posture -1;
%        targetHoldTarget = 1;
%        tempTrajStruct = targetHoldTrajStruct([targetHoldTrajStruct.posture]==targetHoldPosture & [targetHoldTrajStruct.target]==targetHoldTarget);
%        for i = 1:size(tempTrajStruct,2)
%            for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
%                trialAvg = tempTrajStruct(i).allSmoothFR(j).trialAvg;
%                allObsTH = vertcat(allObsTH,trialAvg);
%            end
%        end
%        minNumSamples = min(size(allObsCH,1),size(allObsTH,1));
%         allObs = vertcat(allObs,datasample(allObsCH,minNumSamples,'Replace',false));
%         allObs = vertcat(allObs,datasample(allObsTH,minNumSamples,'Replace',false));
%        obsStruct(structInd).allObs = allObs;
%        obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
%        structInd = structInd + 1;
%     end
%     [postureLDA] = doLDA(obsStruct);
%     
    
%% Find peri-reach target plane 
    % Get minimum number of timestamps in condition averages
    numTimestamps = [];
    for i = 1:size(reachTrajStruct,2)
        numTimestamps(i) = length(reachTrajStruct(i).avgSmoothFR.timestamps);
    end
    figure
    histogram(numTimestamps)
    xlabel('Number of 25ms bins')
    ylabel('Number of trials')
    [minNumTimestamps,i] = min(numTimestamps);
        
    %Form X, containing trial-averaged data for each condition
    X = NaN(minNumTimestamps,numTargets,numPostures,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgSmoothFR.traj;
            X(:,targetInd,postureInd,:) = traj(1:minNumTimestamps,:); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
        
    %Perform marginalizations and store them
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Offsets
    targetMargOffset = mean(X,[1 3],'omitnan') - CIMargOffset;
    postureMargOffset = mean(X,[1 2],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj - targetMargOffset - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - postureMargOffset - CIMargOffset;
    targetTrajNoCP = X - CIMargTraj - CIMargOffset - postureMargTraj - postureMargOffset;
    postureTrajNoC = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
            
    %Do PCA on Marginalizations    
    CIMargTraj = squeeze(CIMargTraj);
        CIMargOffset = squeeze(CIMargOffset);
        [CIPCA,score,latent,tsquared,explained,mu] = pca(CIMargTraj); 
    targetMargTraj = squeeze(targetMargTraj);
        targetMargTraj = reshape(targetMargTraj,[numTargets*minNumTimestamps,numChannels]);
        [targetPCA,score,latent,tsquared,explained,mu] = pca(targetMargTraj); 
    targetTrajNoCP = squeeze(targetTrajNoCP);
        targetTrajNoCP = reshape(targetTrajNoCP,[numTargets*numPostures*minNumTimestamps,numChannels]);
        [targetNoCPPCA,score,latent,tsquared,explained,mu] = pca(targetTrajNoCP); 
    postureMargOffset = squeeze(postureMargOffset);
        [posturePCA,score,latent,tsquared,explained,mu] = pca(postureMargOffset); 
    postureTrajNoC = squeeze(postureTrajNoC);
        postureTrajNoC = reshape(postureTrajNoC,[numPostures*minNumTimestamps,numChannels]);
        [postureNoCPCA,score,latent,tsquared,explained,mu] = pca(postureTrajNoC); 

            
  %Orthonormalize combinations of axes
    [CPTOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1),CIPCA(:,1)]); CPTOrth = CPTOrth(:,1:3);
    [CPOrth,~] = qr([postureLDA(:,1:2),CIPCA(:,1)]); CPOrth = CPOrth(:,1:3);
    [PTOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
%     [PTOrthLDA,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]); PTOrthLDA = PTOrthLDA(:,1:3);
%     [PTOrthPCA,~] = qr([posturePCA(:,1),targetNoCPPCA(:,1:2)]); PTOrthPCA = PTOrthPCA(:,1:3);
%     [CPTAllOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1:2),CIPCA(:,1)]); CPTAllOrth = CPTAllOrth(:,1:4);

        
     %Add Projections to trajStruct
   % totalVar = trace(cov(allTraj));
    for i = 1:size(trajStruct,2)
        %Add average traces to trajStruct
        %trajStruct(i).PCA.traj = (trajStruct(i).avgSmoothFR.traj-allMu)*allPCA;
%         trajStruct(i).CIPCA.traj = (trajStruct(i).avgSmoothFR.traj-CIMargOffset')*CIPCA;
%         trajStruct(i).targetPCA.traj = trajStruct(i).avgSmoothFR.traj*targetPCA;
%         trajStruct(i).posturePCA.traj = trajStruct(i).avgSmoothFR.traj*posturePCA;
%         trajStruct(i).postureLDA.traj = trajStruct(i).avgSmoothFR.traj*postureLDA;
%         trajStruct(i).CPTOrth.traj = trajStruct(i).avgSmoothFR.traj*CPTOrth;
%         trajStruct(i).CPOrth.traj = trajStruct(i).avgSmoothFR.traj*CPOrth;
%         trajStruct(i).PTOrth.traj = trajStruct(i).avgSmoothFR.traj*PTOrth;
           trajStruct(i).PTOrth.traj = trajStruct(i).avgSmoothFR.traj*PTOrth;
        reachTrajStruct(i).PTOrth.traj = reachTrajStruct(i).avgSmoothFR.traj*PTOrth;
        targetHoldTrajStruct(i).PTOrth.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*PTOrth;
        targetHoldPlotTrajStruct(i).PTOrth.traj = targetHoldPlotTrajStruct(i).avgSmoothFR.traj*PTOrth;
        
        reachTrajStruct(i).postureLDA.traj = reachTrajStruct(i).avgSmoothFR.traj*postureLDA;
        targetHoldTrajStruct(i).postureLDA.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*postureLDA;
%         trajStruct(i).PTOrthLDA.traj = trajStruct(i).avgSmoothFR.traj*PTOrthLDA;
%         trajStruct(i).PTOrthPCA.traj = trajStruct(i).avgSmoothFR.traj*PTOrthPCA;
%         trajStruct(i).CPTAllOrth.traj = trajStruct(i).avgSmoothFR.traj*CPTAllOrth;
%         %Get VAF
%         trajStruct(i).PCA.VAF =  100.*(diag(cov(allTraj*allPCA))')./totalVar;
%         trajStruct(i).PTOrth.VAF =  100.*(diag(cov(allTraj*PTOrth))')./totalVar;
%         trajStruct(i).PTOrthLDA.VAF =  100.*(diag(cov(allTraj*PTOrthLDA))')./totalVar;
%         trajStruct(i).PTOrthPCA.VAF =  100.*(diag(cov(allTraj*PTOrthPCA))')./totalVar;
%         trajStruct(i).postureLDA.VAF = 100.*(diag(cov(allTraj*postureLDA))')./totalVar;
    end
    
    for i = 1:size(centerHoldTrajStruct,2)
        centerHoldTrajStruct(i).PTOrth.traj = centerHoldTrajStruct(i).avgSmoothFR.traj*PTOrth;
    end
    
    
%% Create plots
    %Fig 0 - plot reaching dynamics to get axis limits
    tempTrajStruct = reachTrajStruct;
    figure
    hold on
    timePts = 1:5;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
                traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth.traj(timePts,:); 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end

    view([20 10])
    ax = gca;
    xlims = ax.XLim; ylims = ax.YLim; zlims = ax.ZLim;
zlims = [-2,2];

%% Fig 1 - hold states
    tempTrajStruct = centerHoldTrajStruct;
    figure
    hold on
    %timePts = 1:minNumTimestamps;
    postureHoldPts = []; postureInd = 1;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        if any([tempTrajStruct.posture]==posture)
            traj = tempTrajStruct([tempTrajStruct.posture]==posture).PTOrth.traj; 
            traj = mean(traj,1);
            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',40,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            postureHoldPts(postureInd,:) = traj;
            postureInd = postureInd + 1;
        end
    end

    grid on
    view([20 10])
    xlim(xlims); ylim(ylims); zlim(zlims);
    xticklabels({}); yticklabels({}); zticklabels({}); 
    xlabel('Target Dim 1')
    ylabel('Target Dim 2') 
    zlabel('Posture Dim 1')
    
    
%% Fig 2 - traj emanating from hold states
    tempTrajStruct = reachTrajStruct;
    figure
    hold on
    timePts = 1:3;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
                traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth.traj(timePts,:); 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end

    grid on
    view([20 10])
    xlim(xlims); ylim(ylims); zlim(zlims);
    xticklabels({}); yticklabels({}); zticklabels({}); 
    xlabel('Target Dim 1')
    ylabel('Target Dim 2') 
    zlabel('Posture Dim 1')

    
%% Fig 4 - traj emanating from hold states
    tempTrajStruct = trajStruct;
    figure
    hold on
    %timePts = 1:3;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
                traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth.traj; 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end

    grid on
    view([20 10])
    %xlim(xlims); ylim(ylims); zlim(zlims);
    xticklabels({}); yticklabels({}); zticklabels({}); 
    xlabel('Target Dim 1')
    ylabel('Target Dim 2') 
    zlabel('Posture Dim 1')

    %% Fig 3 - traj returning to hold states
    tempTrajStruct = targetHoldPlotTrajStruct;
    figure
    hold on
    xDim = 2; yDim = 3; zDim = 1;
    for posture = [1:4]
        for target = [1,5]
            if posture == 1 && target ==5
            elseif posture == 4 && target == 1
            else
                if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
                    traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth.traj; 

                    %if
                    plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));

                    %Plot posture hold pts
                    traj = postureHoldPts(posture,:);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',40,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                end
            end
        end
    end
    
    grid on
    view([20 10])
    xlim(xlims); ylim(ylims); %zlim(zlims);
    xticklabels({}); yticklabels({}); zticklabels({}); 
    xlabel('Target Dim 1')
    ylabel('Target Dim 2') 
    zlabel('Posture Dim 1')
    
    
    
% %% Fig 3 - traj returning to hold states
%     tempTrajStruct = targetHoldTrajStruct;
%     figure
%     hold on
%     xDim = 1; yDim = 2; zDim = 3;
%     for posture = [1:4]
%         for target = [1,5]
%             if posture == 1 && target ==5
%             elseif posture == 4 && target == 1
%             else
%                 if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
%                     traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth.traj; 
%                     plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
%                     plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
%                     plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
% 
%                     %Plot posture hold pts
%                     traj = postureHoldPts(posture,:);
%                     plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',40,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
%                 end
%             end
%         end
%     end
%     
%     grid on
%     view([20 10])
%     xlim(xlims); ylim(ylims); %zlim(zlims);
%     xticklabels({}); yticklabels({}); zticklabels({}); 
%     xlabel('Posture Dim 1')
%     ylabel('Posture Dim 2') 
%     zlabel('Posture Dim 3')
    
%% Local function for performing LDA
     function [LDAproj] = doLDA(obsStruct)
%         %Preallocate
%         minNumObs = min([obsStruct.numObs]);
        numClasses = size(obsStruct,2);
        numDims = size(obsStruct(1).allObs,2);
% %         obs = NaN(minNumObs*numClasses,numDims); labels =  NaN(minNumObs*numClasses,1); 
% 
%         %Fill
%         k = 1;
%         for i = 1:size(obsStruct,2)
%             totalClassObs = size(obsStruct(i).allObs,1);
%             obs(k:k+minNumObs-1,:) = obsStruct(i).allObs(randsample(totalClassObs,minNumObs),:);
%             labels(k:k+minNumObs-1,1) = ones(minNumObs,1).*obsStruct(i).label;
%             k = k+minNumObs;
%         end
        obs = []; labels = [];
        for i = 1:numel(obsStruct)
            obs = vertcat(obs,obsStruct(i).allObs);
            labels = vertcat(labels,obsStruct(i).label*ones(size(obsStruct(i).allObs,1),1));
        end
        
        LDAproj = fisherLDA(obs, labels);
        [LDAproj,~] = qr(LDAproj);
        LDAproj = LDAproj(:,1:numClasses-1);
     end   