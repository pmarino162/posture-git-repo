clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230201\Fig 5';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = orli;
    tcmap = customRainbow;
    mcmap = tcmap([1,4,6],:);
    
%% Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],'result',[]);
    structInd = 1;
    
%% Parameters
   %manVarThreshold = 90;
   numManPCs = 10;
   numIterations = 1;

%% Run loop for each dataset
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};
    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};
    
    isoDatasetList = {'E20200116','E20200117','E20200120'};
    
    %task = 'bci';
    task = 'iso';
    for datasetList = {'E20200116'}% reachDatasetList%{'E20210707','N20190226','R20200221'}%{'E20200316','N20171215','R20201020'}%{ 'R20200221'}%bciDatasetList%{'E20210706','E20210707','E20210708','E20210709'}%reachDatasetList
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318','E20200319'}
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
                task = 'bci';
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
                task = 'bci';
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
                task = 'bci';
            %Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
                task = 'reach';
            case {'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
                task = 'reach';
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
                task = 'reach';
            %Iso
            case {'E20200116','E20200117','E20200120'}
                trialInclStates(1).trialName = {'IsometricForce_1D'};   
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',50},{'state','Target','first',250}};
                task = 'iso';
                
        end
        
        %Remove any long or short trials - Reach
        if strcmpi(task,'reach')
            kinData = [Data.kinData];
            rxnTime = [kinData.rxnTime];
            reachTime = [kinData.reachTime];
            figure
                histogram(rxnTime)
                xlabel('Reaction Time')
                ylabel('Number of conditions')
            figure
                histogram(reachTime)
                xlabel('Reach Time')
                ylabel('Number of conditions')
           rxnCutoff = prctile(rxnTime,95);
           reachCutoff = prctile(reachTime,95);
           rmTrials =  unique([find(rxnTime > rxnCutoff), find(reachTime > reachCutoff)]);
           Data(rmTrials) = [];
       end

       trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
      
        %Remove long/short trials, BCI
        if strcmpi(task,'bci')
            reachTime = [];
            for i = 1:size(trajStruct,2)
               numTraj = size(trajStruct(i).allSmoothFR,2);
               for j = 1:numTraj
                  reachTime = [reachTime,size(trajStruct(i).allSmoothFR(j).traj,1)]; 
               end
            end     
            figure
                histogram(reachTime)
                xlabel('Reach Time')
                ylabel('Number of trials')
           reachCutoff = prctile(reachTime,95);
           for i = 1:size(trajStruct,2)
               numTraj = size(trajStruct(i).allSmoothFR,2);
               reachTime = [];
               for j = 1:numTraj
                  reachTime = [reachTime,size(trajStruct(i).allSmoothFR(j).traj,1)]; 
               end
               rmTrials =  find(reachTime > reachCutoff);
               trajStruct(i).allSmoothFR(rmTrials) = [];
            end  
        end
        
        
%         %Remove long/short trials, Iso
%         if strcmpi(task,'iso')
%             reachTime = [];
%             for i = 1:size(trajStruct,2)
%                numTraj = size(trajStruct(i).allSmoothFR,2);
%                for j = 1:numTraj
%                   reachTime = [reachTime,size(trajStruct(i).allSmoothFR(j).traj,1)]; 
%                end
%             end     
%             figure
%                 histogram(reachTime)
%                 xlabel('Reach Time')
%                 ylabel('Number of trials')
%            reachCutoff = prctile(reachTime,95);
%            for i = 1:size(trajStruct,2)
%                numTraj = size(trajStruct(i).allSmoothFR,2);
%                reachTime = [];
%                for j = 1:numTraj
%                   reachTime = [reachTime,size(trajStruct(i).allSmoothFR(j).traj,1)]; 
%                end
%                rmTrials =  find(reachTime > reachCutoff);
%                trajStruct(i).allSmoothFR(rmTrials) = [];
%             end  
%         end
 
        % Get number of trials for each condition, numTimestamps in each
        % trail
        numCondTraj = [];
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allSmoothFR,2);
           numCondTraj = [numCondTraj,numTraj];
        end        
        figure
        histogram(numCondTraj)
        xlabel('Number of trials')
        ylabel('Number of conditions')

        %Remove any conditions for which there weren't enough trials
        cutoffNumTraj = 10;
        trajStruct = trajStruct(numCondTraj >= cutoffNumTraj);
        
        %Keep only postures with all targets
        postureList = unique([trajStruct.posture]);
        targetList = unique([trajStruct.target]); 
        keepPosture = [];
        for posture = postureList
            tempTrajStruct = trajStruct([trajStruct.posture]==posture);
            postureTargetList = [tempTrajStruct.target];
            if isequal(postureTargetList,targetList)
                keepPosture = [posture,keepPosture];
            end
        end
        trajStruct = trajStruct(ismember([trajStruct.posture],keepPosture));
        
        %Get minimum number of trials and timestamps
        numTimestamps = [];
        numCondTraj = [];
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allSmoothFR,2);
           numCondTraj = [numCondTraj,numTraj];
           numTimestamps = [numTimestamps,size(trajStruct(i).avgSmoothFR.timestamps,2)]; 
        end
        [minNumTimestamps,i] = min(numTimestamps);
        [minNumCondTraj,i] = min(numCondTraj);
        numPts = minNumTimestamps;
        
        %Manually enter number of points to consider
         switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318','E20200319'}
                numPts = 8;
            case {'N20171215','N20180221'}
                numPts = 12;
            case {'R20201020','R20201021'}
                numPts = 10;
            %Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                numPts = 8;
            case {'R20200221','R20200222'}
                numPts = 8;
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                numPts = 8;
             %Iso
            case {'E20200116','E20200117','E20200120'}
                numPts = 9;
        end
        
        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);

        %Project all data down to top PCs
        allTraj = NaN(numConditions*minNumTimestamps,numChannels);
        j = 1;
        for i = 1:numConditions
           allTraj(j:j+minNumTimestamps-1,:) = trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps,:);
           j = j + minNumTimestamps;
        end
        [allPCs,~,~,~,explained,allMu] = pca(allTraj);
%         numManPCs = 1;
%         while sum(explained(1:numManPCs)) < manVarThreshold
%             numManPCs = numManPCs + 1;
%         end
        for i = 1:size(trajStruct,2)
           trajStruct(i).avgSmoothFR.traj = (trajStruct(i).avgSmoothFR.traj-allMu)*allPCs(:,1:numManPCs); 
           for j = 1:size(trajStruct(i).allSmoothFR,2)
                trajStruct(i).allSmoothFR(j).traj = (trajStruct(i).allSmoothFR(j).traj-allMu)*allPCs(:,1:numManPCs); 
           end
        end
        numDims = numManPCs;
        
        %Set up sessionResultStruct
        sessionResultStruct = struct('predictedPosture',[],'predictedTarget',[],'predictorPosture',[],'predictorTarget',[],... 
            'errorBeforeShift',[],'errorAfterShift',[],'normErrorBeforeShift',[],'normErrorAfterShift',[],...
            'R2BeforeShift',[],'R2AfterShift',[]);
        sessionResultStructInd = 1;
        for predictedPosture = postureList
            for predictedTarget = targetList
                for predictorPosture = postureList
                    sessionResultStruct(sessionResultStructInd).predictedPosture = predictedPosture;
                    sessionResultStruct(sessionResultStructInd).predictedTarget = predictedTarget;
                    sessionResultStruct(sessionResultStructInd).predictorPosture = predictorPosture;
                    sessionResultStruct(sessionResultStructInd).predictorTarget = predictedTarget;
                    sessionResultStruct(sessionResultStructInd).errorBeforeShift = NaN(1,numIterations);
                    sessionResultStruct(sessionResultStructInd).errorAfterShift = NaN(1,numIterations);
                    sessionResultStruct(sessionResultStructInd).normErrorBeforeShift = NaN(1,numIterations);
                    sessionResultStruct(sessionResultStructInd).normErrorAfterShift = NaN(1,numIterations);
                    sessionResultStruct(sessionResultStructInd).R2BeforeShift = NaN(1,numIterations);
                    sessionResultStruct(sessionResultStructInd).R2AfterShift = NaN(1,numIterations);
                    sessionResultStructInd = sessionResultStructInd + 1;
                end
            end
        end
        
        %Set up structs of predicted and actual data for visualization
        %of sample model performance
        actualStruct = struct('target',[],'posture',[],'traj',[]);
        predictedStruct = struct('target',[],'posture',[],'traj',[]);
        structInd = 1;
        for posture = postureList
            for target = targetList
                actualStruct(structInd).posture = posture;
                actualStruct(structInd).target = target;
                predictedStruct(structInd).posture = posture;
                predictedStruct(structInd).target = target;
                structInd = structInd + 1;
            end
        end
        
        numIterations = 1;
        %Do computation
        numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw
        for i = 1:numIterations
            i
            %Split trajStruct into two groups
            trajStruct1 = trajStruct;
            trajStruct2 = trajStruct;
            for j = 1:size(trajStruct,2)
                numTraj = size(trajStruct(j).allSmoothFR,2);
                sampInd1 = randsample(numTraj,numSample);
                trajStruct1(j).allSmoothFR = trajStruct(j).allSmoothFR(sampInd1);
                [trajStruct1(j).avgSmoothFR.traj,trajStruct1(j).avgSmoothFR.timestamps] = getAvgTraj20211210(trajStruct1(j).allSmoothFR,binWidth);               
                numTrajRemaining = numTraj - numSample;
                sampInd2 = randsample(numTrajRemaining,numSample);
                remainingInd = setdiff(1:numTraj,sampInd1);
                sampInd2 = remainingInd(sampInd2);
                trajStruct2(j).allSmoothFR = trajStruct(j).allSmoothFR(sampInd2);
                [trajStruct2(j).avgSmoothFR.traj,trajStruct2(j).avgSmoothFR.timestamps] = getAvgTraj20211210(trajStruct2(j).allSmoothFR,binWidth);   
            end
            
            %Make every comparison
            for predictedPosture = postureList
                for predictedTarget = targetList
                    predictedCond = find([trajStruct1.target]==predictedTarget & [trajStruct1.posture]==predictedPosture);
                    traj1 = trajStruct1(predictedCond).avgSmoothFR.traj;
                    actualStruct([actualStruct.posture]==predictedPosture & [actualStruct.target]==predictedTarget).traj = traj1;
                    %Choose predictor posture to store for visualization in advance
                    storePosture = randsample(setdiff(postureList,predictedPosture),1);
                    for predictorPosture = postureList
                        predictorCond = find([trajStruct2.target]==predictedTarget & [trajStruct2.posture]==predictorPosture);
                        traj2 = trajStruct2(predictorCond).avgSmoothFR.traj;
                        unshiftTraj2 = traj2;
                        if min([size(traj1,1),size(traj2,1)]) < numPts
                            tempNumPts = min([size(traj1,1),size(traj2,1)]);
                        else
                            tempNumPts = numPts;
                        end
                        %cvShift - always 0 for within posture comparisons
                        cvShift = getCVShift(predictedPosture,predictorPosture,predictedTarget,tempNumPts,numDims,numTargets,trajStruct1,trajStruct2,targetList);
                        shiftTraj2 = traj2 + cvShift;
                        
                        %Store predicted traj for visualization
                        if predictorPosture == storePosture
                            predictedStruct([predictedStruct.posture]==predictedPosture & [predictedStruct.target]==predictedTarget).traj = shiftTraj2;
                        end
                        
                        %Compute absolute errors
                        errorBeforeShift = getMeanDist(traj1,unshiftTraj2,tempNumPts);
                        errorAfterShift = getMeanDist(traj1,shiftTraj2,tempNumPts);
                        R2BeforeShift = getR2(traj1,unshiftTraj2,tempNumPts);
                        R2AfterShift = getR2(traj1,shiftTraj2,tempNumPts);
                        %Store results in sessionResultStruct; compute normalized errors
                        sessionResultStructInd = find([sessionResultStruct.predictedPosture]==predictedPosture & ...
                            [sessionResultStruct.predictedTarget]==predictedTarget & [sessionResultStruct.predictorPosture]==predictorPosture);
                        sessionResultStruct(sessionResultStructInd).errorBeforeShift(i) = errorBeforeShift;
                        sessionResultStruct(sessionResultStructInd).errorAfterShift(i) = errorAfterShift;                       
                       
                        sessionResultStruct(sessionResultStructInd).R2BeforeShift(i) = R2BeforeShift;
                        sessionResultStruct(sessionResultStructInd).R2AfterShift(i) = R2AfterShift;
                    end
                end
            end
            
            %Get normalized errors
            for predictedPosture = postureList
                for predictedTarget = targetList
                    for predictorPosture = postureList
                        sessionResultStructInd = find([sessionResultStruct.predictedPosture]==predictedPosture & ...
                            [sessionResultStruct.predictedTarget]==predictedTarget & [sessionResultStruct.predictorPosture]==predictorPosture);
                        baselineError = sessionResultStruct([sessionResultStruct.predictedPosture]==predictedPosture & ...
                            [sessionResultStruct.predictedTarget]==predictedTarget & [sessionResultStruct.predictorPosture]==predictedPosture).errorAfterShift(i);
                        errorBeforeShift = sessionResultStruct(sessionResultStructInd).errorBeforeShift(i);
                        errorAfterShift = sessionResultStruct(sessionResultStructInd).errorAfterShift(i);
                        normErrorBeforeShift = errorBeforeShift./baselineError;
                        normErrorAfterShift = errorAfterShift./baselineError;
                        sessionResultStruct(sessionResultStructInd).normErrorBeforeShift(i) = normErrorBeforeShift;
                        sessionResultStruct(sessionResultStructInd).normErrorAfterShift(i) = normErrorAfterShift;
                    end
                end
            end
            
            
        end

        %Add mean
        for sessionResultStructInd = 1:numel(sessionResultStruct)
           sessionResultStruct(sessionResultStructInd).meanNormErrorAfterShift = mean([sessionResultStruct(sessionResultStructInd).normErrorAfterShift]);
        end
        
        %Remove self-comparisons & repeat comparisons
        rmInd = [];
        for sessionResultStructInd = 1:numel(sessionResultStruct)
           predictedPosture = sessionResultStruct(sessionResultStructInd).predictedPosture;
           predictorPosture = sessionResultStruct(sessionResultStructInd).predictorPosture;
           if predictedPosture==predictorPosture
               rmInd = [sessionResultStructInd,rmInd];
           end
           if predictedPosture > predictedPosture
              rmInd = [sessionResultStructInd,rmInd];
           end
        end
        sessionResultStruct(rmInd) = [];
        
        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).result = sessionResultStruct;
        structInd = structInd + 1;
    end
    
%% Plot results
%Get minTimestamps from any condition (predicted or actual)
    % Get timestamps dist and min number timestamps
    numTimestamps = [];
    for i = 1:size(actualStruct,2)
        numTimestamps = [numTimestamps,size(actualStruct(i).traj,1),size(predictedStruct(i).traj,1)]; 
    end
    [minNumTimestamps,~] = min(numTimestamps);
    if minNumTimestamps < numPts
        numPts = minNumTimestamps;
    end
    
%Fit dPCA to actual data
    %Form X, containing trial-averaged data for each condition
    X = NaN(numPts,numTargets,numPostures,numDims);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = actualStruct([actualStruct.posture]==posture & [actualStruct.target]==target).traj;
            X(:,targetInd,postureInd,:) = traj(1:numPts,:); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end

    Xdpca = permute(X,[4,3,2,1]);
    combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
    margNames = {'P', 'T', 'CI', 'PTI'};        

    [W,V,whichMarg] = dpca(Xdpca, numManPCs, ...
        'combinedParams', combinedParams);
    explVar = dpca_explainedVariance(Xdpca, W, V, ...
        'combinedParams', combinedParams);   

    postureMargID = find(strcmp(margNames,'P'));
    targetMargID = find(strcmp(margNames,'T'));
    postureDims = find(whichMarg==postureMargID,2);
    targetDims = find(whichMarg==targetMargID,2);

    %Orthonormalize pdpca1 and tdpca1-2
    pDPCA = W(:,postureDims);
    tDPCA = W(:,targetDims);
    [PTOrth,~] = qr([pDPCA(:,1),tDPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
        
%Plot actual data, save plot lims
    figure; hold on; fs = 14;
    if numel(postureList)>5
        postureList = postureList(1:5);
    end

    trajStruct = actualStruct; 

    xDim = 2; yDim = 3; zDim = 1;
    postureInd = 1;
    for posture = postureList
        for target = targetList
            if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj(1:numPts,:)*PTOrth; 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(postureInd,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(postureInd,:),'MarkerFaceColor',pcmap(postureInd,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(postureInd,:),'MarkerFaceColor',pcmap(postureInd,:));
            end
        end
        postureInd = postureInd +1;
    end       


    xlabel('Goal Dim 1'); ylabel('Goal Dim 2') ;zlabel('Posture Dim 1')
    view([20 10])
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    ax = gca;
    xlimits = ax.XLim; ylimits = ax.YLim; zlimits = ax.ZLim; 
    
    xlimits = [-6,xlimits(2)];
    xlim(xlimits);
    
    xticksPlot = ax.XTick; yticksPlot = ax.YTick; zticksPlot = ax.ZTick; 
    xticklabels({}); yticklabels({}); zticklabels({}); 
    grid on

    if saveFig
        saveas(gcf,fullfile(saveDir,[task,'_actualExample.svg']));
    end
        
%Plot predicted data 
    figure;hold on; fs = 14;
        trajStruct = predictedStruct; 
        
        xDim = 2; yDim = 3; zDim = 1;
        postureInd = 1;
        for posture = postureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj(1:numPts,:)*PTOrth; 
                    plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(postureInd,:),'LineWidth',2);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(postureInd,:),'MarkerFaceColor',pcmap(postureInd,:));
                    plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(postureInd,:),'MarkerFaceColor',pcmap(postureInd,:));
                end
            end
            postureInd = postureInd +1;
        end       
        
        

        xlabel('Goal Dim 1'); ylabel('Goal Dim 2') ;zlabel('Posture Dim 1')
        view([20 10])
        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        ax = gca;
        xlim(xlimits); ylim(ylimits); zlim(zlimits);
        xticks(xticksPlot); yticks(yticksPlot); zticks(zticksPlot);
        xticklabels({}); yticklabels({}); zticklabels({}); 
        grid on
             
        if saveFig
            saveas(gcf,fullfile(saveDir,[task,'_modelPredExample.svg']));
        end
      
       


        

        


%% Local function   
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
    
    %Get R2
    function r2 = getR2(traj1,traj2,numPts)
        SST = 0;
        SSR = 0;
        traj1mean = mean(traj1(1:numPts,:),1);
        numDims = size(traj1,2);
        for i = 1:numPts
            for j = 1:numDims
                SST = (traj1(i,j)-traj1mean(j)).^2 + SST;
                SSR = (traj1(i,j)-traj2(i,j)).^2 + SSR;
            end
        end
        r2 = 1-SSR/SST;
    end
    
    %Get cvShift
    function cvShift = getCVShift(posture1,posture2,curTarget,numPts,numDims,numTargets,trajStruct1,trajStruct2,targetList)  
        %Get optimal alpha for every target besides current one
        alphaList = NaN(numTargets-1,numDims);
        listInd = 1;
        for target = setdiff(targetList,curTarget)
            traj1 = trajStruct2([trajStruct2.posture]==posture1 & [trajStruct2.target]==target).avgSmoothFR.traj;
            traj2 = trajStruct2([trajStruct2.posture]==posture2 & [trajStruct2.target]==target).avgSmoothFR.traj;
            if min([size(traj1,1),size(traj2,1)]) < numPts
                tempNumPts = min([size(traj1,1),size(traj2,1)]);
            else
                tempNumPts = numPts;
            end        
            alphaList(listInd,:) = getOptimalAlpha(traj1,traj2,tempNumPts);
            listInd = listInd + 1;
        end
        %Average together
        cvShift = mean(alphaList,1);
    end
    
   
    
