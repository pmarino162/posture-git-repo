clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'D:\Posture Paper\20221123\Figures\Figure 5 - similar traj';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = orli;
    scmap = customRainbow;
    APcolor = scmap(4,:);
    ATcolor = scmap(8,:);
    refColor = [0.5 0.5 0.5];
    fs = 14;
    
%% Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],'vPP',[],'vPT',[],'vTT',[],'vTP',[],'APT',[],'ATP',[],'nullAPT',[],'nullATP',[]);
    structInd = 1;

%% Parameters
    % How many PC's to describe the 'manifold'? 
    manVarThreshold = 90;
    
%% Get traj struct
    bciDatasetList = {'E20200316','E20200317','E20200318','N20171215','N20180221','R20201020'};
    task = 'BCI';

    datasetList = {'E20200116'};
    dataset = datasetList{1,1};
    [Data,zScoreParams] = loadData(dataset);
        
    %Get trajStruct
    binWidth = 25;
    kernelStdDev = 25;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    switch dataset
        %BCI
        case {'E20200316','E20200317','E20200318'}
            trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
            condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}}; 
        case {'N20171215','N20180221'}
            trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
        case {'R20201020','R20201021'}
            trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
        case 'E20210901'
            taskID = [Data.conditionData]; taskID = [taskID.taskID];
            Data = Data(taskID==1);
            trialInclStates(1).trialName = {'BCI Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
        %Iso
        case {'E20200116','E20200117','E20200120'}
            trialInclStates(1).trialName = {'IsometricForce_1D'};   
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',200}};

        %Reaching
        case{'E20210706','E20210707','E20210708','E20210709','E20210710'}
            trialInclStates(1).trialName = {'GridReaching'};  
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
        case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
            trialInclStates(1).trialName = {'Nigel Dissociation'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
        case{'R20200221','R20200222'}
            trialInclStates(1).trialName = {'Rocky Dissociation'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};

    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

% %% Eliminate trials w unusual length
%     numTimestamps = [];
%     for i = 1:size(trajStruct,2)
%        numTraj = size(trajStruct(i).allSmoothFR,2);
%        for j = 1:numTraj
%           numTimestamps = [numTimestamps,size(trajStruct(i).allSmoothFR(j).timestamps,2)]; 
%        end
%     end
%     histogram(numTimestamps,[1:max(numTimestamps)])
%     %Eliminate trials outside of range
%     lowerBound = prctile(numTimestamps,5);
%     upperBound = prctile(numTimestamps,90);
    
%% Get timestamps dist and min number timestamps
    numTimestamps = [];
    numCondTraj = [];
    for i = 1:size(trajStruct,2)
       numTraj = size(trajStruct(i).allSmoothFR,2);
       numCondTraj = [numCondTraj,numTraj];
       numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    

    [minNumTimestamps,i] = min(numTimestamps);
    [minNumCondTraj,i] = min(numCondTraj);
    
    %Get posture & target lists
    postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2); numConditions = size(trajStruct,2);

    %Keep only postures with all targets for earl
    if strcmpi(dataset(1),'E')
        keepPosture = [];
        for posture = postureList
            tempTrajStruct = trajStruct([trajStruct.posture]==posture);
            postureTargetList = [tempTrajStruct.target];
            if isequal(postureTargetList,targetList)
                keepPosture = [posture,keepPosture];
            end
        end
        trajStruct = trajStruct(ismember([trajStruct.posture],keepPosture));
    end

    %Update posture list
    postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
 
%% Get CPT space for visualization
    %Form X
    X = NaN(minNumTimestamps,numTargets,numPostures,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            X(:,targetInd,postureInd,:) = traj(1:minNumTimestamps,:); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
            
    %Compute within space variances
    allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
    [allPCs,~,~,~,explained,allMu] = pca(allTraj);
            
    % Do marginalizations of X
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;

    %Compute PCs of marginalizations
    CIMargTraj = squeeze(CIMargTraj);
        CIMargOffset = squeeze(CIMargOffset);
        [CIPCs,~,~,~,explainedCI,CIMu] = pca(CIMargTraj); 
    targetMargTraj = squeeze(targetMargTraj);
        targetMargTraj = reshape(targetMargTraj,[numTargets*minNumTimestamps,numChannels]);
        [targetPCs,~,~,~,explainedT,targetMargMu] = pca(targetMargTraj); 
    postureMargTraj = squeeze(postureMargTraj);
        postureMargTraj = reshape(postureMargTraj,[numPostures*minNumTimestamps,numChannels]);
        [posturePCs,~,~,~,explainedP,postureMargMu] = pca(postureMargTraj); 

    %Form joint subspace for full traj visualization      
    [CPTOrth,~] = qr([CIPCs(:,1),posturePCs(:,1),targetPCs(:,1)]); CPTOrth = CPTOrth(:,1:3);  
    
%% Do heavy lifting
  % Create distance distributions
    resultStruct = struct('cond1P',[],'cond1T',[],'cond2P',[],'cond2T',[],'dist',[],'alpha',[],'fullAlpha',[],'uDist',[],'uAlpha',[],'uFullAlpha',[],'traj1Sample',[],'traj2Sample',[],'alphaSample',[],'distSample',[]);
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
            saveSampleNumber = randi(numDraws); %Draw number to save for visualization
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
                numPts = min([size(traj1,1),size(traj2,1)]);
                alpha = getOptimalAlpha(traj1,traj2,numPts);
%                 alpha = traj1(1,:)-traj2(1,:);
                %Shift traj2 by alpha
                unshiftTraj2 = traj2;
                traj2 = traj2 + alpha;
                resultStruct(structInd).fullAlpha(i,:) = alpha;
                fullAlpha = alpha;
                alpha = vecnorm(alpha);
                resultStruct(structInd).alpha(i) = alpha;
                %Get residual distance
                dist = getMeanDist(traj1,traj2,numPts);
                resultStruct(structInd).dist(i) = dist;
                %Save sample 
                if i==saveSampleNumber
                    resultStruct(structInd).traj1Sample = traj1;
                    resultStruct(structInd).traj2Sample = unshiftTraj2;
                    resultStruct(structInd).shiftTraj2Sample = traj2;
                    resultStruct(structInd).alphaSample = alpha;
                    resultStruct(structInd).distSample = dist;           
                end
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
    %Across posture comparisons
    acrossPostureResult = resultStruct([cond1P ~= cond2P] & [cond1T == cond2T]);
    acrossPostureDist = [withinConditionDist,[acrossPostureResult.uDist]];
    acrossPostureAlpha = [withinConditionAlpha,[acrossPostureResult.uAlpha]];
    rPosture = corrcoef([acrossPostureAlpha',acrossPostureDist']);
    rPosture = rPosture(1,2);
    XPosture = [ones(length(acrossPostureAlpha),1),acrossPostureAlpha'];
    YPosture = acrossPostureDist';
    BPosture = inv(XPosture'*XPosture)*XPosture'*YPosture;

    %Across target comparisons
    acrossTargetResult = resultStruct([cond1P == cond2P] & [cond1T ~= cond2T]);
    acrossTargetDist = [withinConditionDist,[acrossTargetResult.uDist]];
    acrossTargetAlpha = [withinConditionAlpha,[acrossTargetResult.uAlpha]];
    rTarget = corrcoef([acrossTargetAlpha',acrossTargetDist']);
    rTarget = rTarget(1,2);
    XTarget = [ones(length(acrossTargetAlpha),1),acrossTargetAlpha'];
    YTarget = acrossTargetDist';
    BTarget = inv(XTarget'*XTarget)*XTarget'*YTarget;

    f = figure; f.Position = [40 40 1000 750];
    hold on
    for i = 1:length(acrossTargetAlpha)
        plot(acrossTargetAlpha(i),acrossTargetDist(i),'o','MarkerSize',10,...
            'MarkerEdgeColor',ATcolor,'MarkerFaceColor',ATcolor);
    end
    for i = 1:length(acrossPostureAlpha)
           plot(acrossPostureAlpha(i),acrossPostureDist(i),'o','MarkerSize',10,...
            'MarkerEdgeColor',APcolor,'MarkerFaceColor',APcolor);

    end
    ax=gca;
    xLim = ax.XLim;
    temp = [ones(2,1),xLim'];
    plot(xLim,temp*BPosture,'Color',APcolor,'LineWidth',1)
    plot(xLim,temp*BTarget,'Color',ATcolor,'LineWidth',1)
     set(gca, 'TickDir', 'out')
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    xlabel('Applied Shift (SD)')
    ylabel('Remaining Trajectory Difference (SD)')
    yLim = ax.YLim;
    yRange = yLim(2)-yLim(1); 

  
    if saveFig
       saveas(gcf,fullfile(saveDir,[dataset,'_Scatter.svg']));
    end
    
%% Visualize random samples
    %Across posture
        cond1P = 1; cond1T = 1;
        cond2P = 2; cond2T = 1;
        structInd = find([resultStruct.cond1P] == cond1P & [resultStruct.cond1T] == cond1T & [resultStruct.cond2P] == cond2P & [resultStruct.cond2T] == cond2T);
        traj1 = resultStruct(structInd).traj1Sample(1:minNumTimestamps,:);
        traj2 = resultStruct(structInd).traj2Sample(1:minNumTimestamps,:);
        traj2Shifted = resultStruct(structInd).shiftTraj2Sample(1:minNumTimestamps,:);
            
        traj1ProjAP = traj1*CPTOrth;
        traj2ProjAP = traj2*CPTOrth;
        traj2ShiftedProjAP = traj2Shifted*CPTOrth;
        numPts = min([size(traj1,1),size(traj2,1)]);
        distBeforeAP = getMeanDist(traj1,traj2,numPts);
        alphaSampleAP = resultStruct(structInd).alphaSample;
        distAfterAP = resultStruct(structInd).distSample;
    
    %Across target
        cond1P = 1; cond1T = 1;
        cond2P = 1; cond2T = 5;
        structInd = find([resultStruct.cond1P] == cond1P & [resultStruct.cond1T] == cond1T & [resultStruct.cond2P] == cond2P & [resultStruct.cond2T] == cond2T);
        traj1 = resultStruct(structInd).traj1Sample(1:minNumTimestamps,:);
        traj2 = resultStruct(structInd).traj2Sample(1:minNumTimestamps,:);
        traj2Shifted = resultStruct(structInd).shiftTraj2Sample(1:minNumTimestamps,:);
        
        traj1ProjAT = traj1*CPTOrth;
        traj2ProjAT = traj2*CPTOrth;
        traj2ShiftedProjAT = traj2Shifted*CPTOrth;
        numPts = min([size(traj1,1),size(traj2,1)]);
        distBeforeAT = getMeanDist(traj1,traj2,numPts);
        alphaSampleAT = resultStruct(structInd).alphaSample;
        distAfterAT = resultStruct(structInd).distSample;
                        
    %Get axis limits
        figure; hold on
        plot3(traj1ProjAP(:,1),traj1ProjAP(:,2),traj1ProjAP(:,3),'LineWidth',2,'Color','r');
        plot3(traj2ProjAP(:,1),traj2ProjAP(:,2),traj2ProjAP(:,3),'LineWidth',2,'Color','r');
        plot3(traj2ShiftedProjAP(:,1),traj2ShiftedProjAP(:,2),traj2ShiftedProjAP(:,3),'LineWidth',2,'Color','r');
        plot3(traj1ProjAT(:,1),traj1ProjAT(:,2),traj1ProjAT(:,3),'LineWidth',2,'Color','r');
        plot3(traj2ProjAT(:,1),traj2ProjAT(:,2),traj2ProjAT(:,3),'LineWidth',2,'Color','r');
        plot3(traj2ShiftedProjAT(:,1),traj2ShiftedProjAT(:,2),traj2ShiftedProjAT(:,3),'LineWidth',2,'Color','r');
        ax = gca;
        xlimits = ax.XLim; ylimits = ax.YLim; zlimits = ax.ZLim;
    
    
        %Unshifted CPT
        figure; hold on
        traj1Proj = traj1ProjAP;
        traj2Proj = traj2ProjAP;
        
        plot3(traj1Proj(:,1),traj1Proj(:,2),traj1Proj(:,3),'LineWidth',2,'Color',refColor);
            plot3(traj1Proj(1,1),traj1Proj(1,2),traj1Proj(1,3),'.','MarkerSize',20,'Color',refColor,'MarkerFaceColor',refColor);
            plot3(traj1Proj(end,1),traj1Proj(end,2),traj1Proj(end,3),'<','MarkerSize',5,'Color',refColor,'MarkerFaceColor',refColor);
        plot3(traj2Proj(:,1),traj2Proj(:,2),traj2Proj(:,3),'LineWidth',2,'Color',APcolor);
            plot3(traj2Proj(1,1),traj2Proj(1,2),traj2Proj(1,3),'.','MarkerSize',20,'Color',APcolor,'MarkerFaceColor',APcolor);
            plot3(traj2Proj(end,1),traj2Proj(end,2),traj2Proj(end,3),'<','MarkerSize',5,'Color',APcolor,'MarkerFaceColor',APcolor);
              
        title(['P',num2str(cond1P),'T',num2str(cond1T),' vs ','P',num2str(cond2P),'T',num2str(cond2T)...
            newline,'Dist = ',num2str(distBeforeAP),newline,'Alpha = ',num2str(alphaSampleAP)])
        ax = gca; ax.XLim = xlimits; ax.YLim = ylimits; ax.ZLim = zlimits;
         axis equal
        grid on;
            
        %Shifted CPT
        traj1Proj = traj1ProjAP;
        traj2Proj = traj2ShiftedProjAP;
        
        figure; hold on
        plot3(traj1Proj(:,1),traj1Proj(:,2),traj1Proj(:,3),'LineWidth',2,'Color',refColor);
            plot3(traj1Proj(1,1),traj1Proj(1,2),traj1Proj(1,3),'.','MarkerSize',20,'Color',refColor,'MarkerFaceColor',refColor);
            plot3(traj1Proj(end,1),traj1Proj(end,2),traj1Proj(end,3),'<','MarkerSize',5,'Color',refColor,'MarkerFaceColor',refColor);
        plot3(traj2Proj(:,1),traj2Proj(:,2),traj2Proj(:,3),'LineWidth',2,'Color',APcolor);
            plot3(traj2Proj(1,1),traj2Proj(1,2),traj2Proj(1,3),'.','MarkerSize',20,'Color',APcolor,'MarkerFaceColor',APcolor);
            plot3(traj2Proj(end,1),traj2Proj(end,2),traj2Proj(end,3),'<','MarkerSize',5,'Color',APcolor,'MarkerFaceColor',APcolor);
        
        title(['P',num2str(cond1P),'T',num2str(cond1T),' vs ','P',num2str(cond2P),'T',num2str(cond2T)...
            newline,'Dist = ',num2str(distAfterAP),newline,'Alpha = ',num2str(alphaSampleAP)])
        ax = gca; ax.XLim = xlimits; ax.YLim = ylimits; ax.ZLim = zlimits;
        axis equal
        grid on;

        %Unshifted CPT
        traj1Proj = traj1ProjAT;
        traj2Proj = traj2ProjAT;
        figure; hold on
        plot3(traj1Proj(:,1),traj1Proj(:,2),traj1Proj(:,3),'LineWidth',2,'Color',refColor);
            plot3(traj1Proj(1,1),traj1Proj(1,2),traj1Proj(1,3),'.','MarkerSize',20,'Color',refColor,'MarkerFaceColor',refColor);
            plot3(traj1Proj(end,1),traj1Proj(end,2),traj1Proj(end,3),'<','MarkerSize',5,'Color',refColor,'MarkerFaceColor',refColor);
        plot3(traj2Proj(:,1),traj2Proj(:,2),traj2Proj(:,3),'LineWidth',2,'Color',ATcolor);
            plot3(traj2Proj(1,1),traj2Proj(1,2),traj2Proj(1,3),'.','MarkerSize',20,'Color',ATcolor,'MarkerFaceColor',ATcolor);
            plot3(traj2Proj(end,1),traj2Proj(end,2),traj2Proj(end,3),'<','MarkerSize',5,'Color',ATcolor,'MarkerFaceColor',ATcolor);
       
        title(['P',num2str(cond1P),'T',num2str(cond1T),' vs ','P',num2str(cond2P),'T',num2str(cond2T)...
            newline,'Dist = ',num2str(distBeforeAT),newline,'Alpha = ',num2str(alphaSampleAT)])
        ax = gca; ax.XLim = xlimits; ax.YLim = ylimits; ax.ZLim = zlimits;
       axis equal
        grid on;
            
        %Shifted CPT
        traj1Proj = traj1ProjAT;
        traj2Proj = traj2ShiftedProjAT;
        figure; hold on
        plot3(traj1Proj(:,1),traj1Proj(:,2),traj1Proj(:,3),'LineWidth',2,'Color',refColor);
            plot3(traj1Proj(1,1),traj1Proj(1,2),traj1Proj(1,3),'.','MarkerSize',20,'Color',refColor,'MarkerFaceColor',refColor);
            plot3(traj1Proj(end,1),traj1Proj(end,2),traj1Proj(end,3),'<','MarkerSize',5,'Color',refColor,'MarkerFaceColor',refColor);
        plot3(traj2Proj(:,1),traj2Proj(:,2),traj2Proj(:,3),'LineWidth',2,'Color',ATcolor);
            plot3(traj2Proj(1,1),traj2Proj(1,2),traj2Proj(1,3),'.','MarkerSize',20,'Color',ATcolor,'MarkerFaceColor',ATcolor);
            plot3(traj2Proj(end,1),traj2Proj(end,2),traj2Proj(end,3),'<','MarkerSize',5,'Color',ATcolor,'MarkerFaceColor',ATcolor);
        
        title(['P',num2str(cond1P),'T',num2str(cond1T),' vs ','P',num2str(cond2P),'T',num2str(cond2T)...
            newline,'Dist = ',num2str(distAfterAT),newline,'Alpha = ',num2str(alphaSampleAT)])
        ax = gca; ax.XLim = xlimits; ax.YLim = ylimits; ax.ZLim = zlimits;
        axis equal
        grid on;


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
    
    
    
    %         %Unshifted PCA
%             traj1Proj = traj1*allPCs;
%             traj2Proj = traj2*allPCs;
%             figure; hold on
%             plot3(traj1Proj(:,1),traj1Proj(:,2),traj1Proj(:,3),'LineWidth',2,'Color','r');
%                 plot3(traj1Proj(1,1),traj1Proj(1,2),traj1Proj(1,3),'.','MarkerSize',20,'Color','r','MarkerFaceColor','r');
%                 plot3(traj1Proj(end,1),traj1Proj(end,2),traj1Proj(end,3),'<','MarkerSize',5,'Color','r','MarkerFaceColor','r');
%             plot3(traj2Proj(:,1),traj2Proj(:,2),traj2Proj(:,3),'LineWidth',2,'Color','b');
%                 plot3(traj2Proj(1,1),traj2Proj(1,2),traj2Proj(1,3),'.','MarkerSize',20,'Color','b','MarkerFaceColor','b');
%                 plot3(traj2Proj(end,1),traj2Proj(end,2),traj2Proj(end,3),'<','MarkerSize',5,'Color','b','MarkerFaceColor','b');
% 
%         %Shifted PCA
%             traj1Proj = traj1*allPCs;
%             traj2Proj = traj2Shifted*allPCs;
%             figure; hold on
%             plot3(traj1Proj(:,1),traj1Proj(:,2),traj1Proj(:,3),'LineWidth',2,'Color','r');
%                 plot3(traj1Proj(1,1),traj1Proj(1,2),traj1Proj(1,3),'.','MarkerSize',20,'Color','r','MarkerFaceColor','r');
%                 plot3(traj1Proj(end,1),traj1Proj(end,2),traj1Proj(end,3),'<','MarkerSize',5,'Color','r','MarkerFaceColor','r');
%             plot3(traj2Proj(:,1),traj2Proj(:,2),traj2Proj(:,3),'LineWidth',2,'Color','b');
%                 plot3(traj2Proj(1,1),traj2Proj(1,2),traj2Proj(1,3),'.','MarkerSize',20,'Color','b','MarkerFaceColor','b');
%                 plot3(traj2Proj(end,1),traj2Proj(end,2),traj2Proj(end,3),'<','MarkerSize',5,'Color','b','MarkerFaceColor','b');
%         