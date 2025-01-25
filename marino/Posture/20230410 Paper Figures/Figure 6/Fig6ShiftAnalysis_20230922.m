clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Error Reduction Figures';
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
   numIterations = 20;

%% Run loop for each dataset
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};
    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};
    
    isoDatasetList = {'E20200116','E20200117','E20200120'};
    
    for datasetList = {'E20200316'}%{'E20200316'}%bciDatasetList% reachDatasetList%{'E20210707','N20190226','R20200221'}%{'E20200316','N20171215','R20201020'}%{ 'R20200221'}%bciDatasetList%{'E20210706','E20210707','E20210708','E20210709'}%reachDatasetList
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
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
                task = 'reach';
            case {'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
                task = 'reach';
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
                task = 'reach';
                
            %Iso
            case {'E20200116','E20200117','E20200120'}
                trialInclStates(1).trialName = {'IsometricForce_1D'};   
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
                task = 'iso';
        end
                
       trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
      
 
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
            'R2BeforeShift',[],'R2AfterShift',[],'pctReduction',[]);
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
                    sessionResultStruct(sessionResultStructInd).pctReduction = NaN(1,numIterations);
                    sessionResultStructInd = sessionResultStructInd + 1;
                end
            end
        end
        
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
                    for predictorPosture = postureList
                        predictedCond = find([trajStruct1.target]==predictedTarget & [trajStruct1.posture]==predictedPosture);
                        predictorCond = find([trajStruct2.target]==predictedTarget & [trajStruct2.posture]==predictorPosture);
                        traj1 = trajStruct1(predictedCond).avgSmoothFR.traj;
                        traj2 = trajStruct2(predictorCond).avgSmoothFR.traj;
                        unshiftTraj2 = traj2;
                        if min([size(traj1,1),size(traj2,1)]) < numPts
                            tempNumPts = min([size(traj1,1),size(traj2,1)]);
                        else
                            tempNumPts = numPts;
                        end
                        %cvShift - always 0 for within posture comparisons
                        alpha = getOptimalAlpha(traj1,traj2,tempNumPts);
%                         if predictedPosture == predictorPosture
%                             shiftTraj2 = traj2; 
%                         else
%                             shiftTraj2 = traj2 + alpha; 
%                         end
                        shiftTraj2 = traj2 + alpha;                        
                        %Compute absolute errors
                        errorBeforeShift = getMeanDist(traj1,unshiftTraj2,tempNumPts);
                        errorAfterShift = getMeanDist(traj1,shiftTraj2,tempNumPts);
                        R2BeforeShift = getR2(traj1,unshiftTraj2,tempNumPts);
                        R2AfterShift = getR2(traj1,shiftTraj2,tempNumPts);
                        %Store results in sessionResultStruct
                        sessionResultStructInd = find([sessionResultStruct.predictedPosture]==predictedPosture & ...
                            [sessionResultStruct.predictedTarget]==predictedTarget & [sessionResultStruct.predictorPosture]==predictorPosture);
                        sessionResultStruct(sessionResultStructInd).errorBeforeShift(i) = errorBeforeShift;
                        sessionResultStruct(sessionResultStructInd).errorAfterShift(i) = errorAfterShift;                                              
                        sessionResultStruct(sessionResultStructInd).R2BeforeShift(i) = R2BeforeShift;
                        sessionResultStruct(sessionResultStructInd).R2AfterShift(i) = R2AfterShift;
                        sessionResultStruct(sessionResultStructInd).pctReduction(i) = 100*(errorBeforeShift-errorAfterShift)/errorBeforeShift;
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
           sessionResultStruct(sessionResultStructInd).meanR2BeforeShift = mean([sessionResultStruct(sessionResultStructInd).R2BeforeShift]);
           sessionResultStruct(sessionResultStructInd).meanR2AfterShift = mean([sessionResultStruct(sessionResultStructInd).R2AfterShift]);
           
           sessionResultStruct(sessionResultStructInd).meanErrorBeforeShift = mean([sessionResultStruct(sessionResultStructInd).errorBeforeShift]);
           sessionResultStruct(sessionResultStructInd).meanErrorAfterShift = mean([sessionResultStruct(sessionResultStructInd).errorAfterShift]);
           
           sessionResultStruct(sessionResultStructInd).meanNormErrorBeforeShift = mean([sessionResultStruct(sessionResultStructInd).normErrorBeforeShift]);
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
           if predictedPosture > predictorPosture
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
    %Across datasets
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
    %Plotting parameters
    fs = 12;
    alpha = 0.3;
    offset = 1/7;   
    sessionOffset = 1/40;
    

    
    %Scatter plot of normalized error
    allMonkeyComparisonStruct = struct('monkey',[],'allMonkeyComparisons',[]);
    f=figure; hold on;
    monkeyInd = 1;
    for monkey = monkeyList
       allMonkeyComparisons = NaN(1,2);
       allMonkeyComparisonInd = 1;
       tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
       for session = 1:numel(tempResultStruct)
           result = tempResultStruct(session).result;
           %Plot normalized error
           for comparison = 1:numel(result)
               %Get before shift and after shift errors
               normErrorBeforeShift = horzcat(result(comparison).normErrorBeforeShift);
               normErrorBeforeShiftMean = mean(normErrorBeforeShift);
               normErrorAfterShift = horzcat(result(comparison).normErrorAfterShift);
               normErrorAfterShiftMean = mean(normErrorAfterShift);
               %Add to allMonkeyComparisons
               allMonkeyComparisons(allMonkeyComparisonInd,:) = [normErrorBeforeShiftMean,normErrorAfterShiftMean];
               allMonkeyComparisonInd = allMonkeyComparisonInd + 1;
               %Plot scatter
               switch monkeyInd
                   case 1
                       scatter(normErrorBeforeShiftMean,normErrorAfterShiftMean,10,'o',...
                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],...
                           'MarkerFaceAlpha',0,'MarkerEdgeAlpha',.5);    
                   case 2
                       scatter(normErrorBeforeShiftMean,normErrorAfterShiftMean,15,'d',...
                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],...
                           'MarkerFaceAlpha',0,'MarkerEdgeAlpha',.5);  
                   case 3
                       scatter(normErrorBeforeShiftMean,normErrorAfterShiftMean,20,'x',...
                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],...
                           'MarkerFaceAlpha',0,'MarkerEdgeAlpha',.5);  
               end
           end
       end
       allMonkeyComparisonsStruct(monkeyInd).monkey = monkey;
       allMonkeyComparisonsStruct(monkeyInd).allMonkeyComparisons = allMonkeyComparisons;
       monkeyInd = monkeyInd + 1;
    end
    monkeyInd = 1;
    for monkey = monkeyList
       allMonkeyComparisons = allMonkeyComparisonsStruct(monkeyInd).allMonkeyComparisons;
       switch monkeyInd
           case 1
                scatter(mean(allMonkeyComparisons(:,1)),mean(allMonkeyComparisons(:,2)),60,'o',...
               'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],...
               'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);    
           case 2           
                   scatter(mean(allMonkeyComparisons(:,1)),mean(allMonkeyComparisons(:,2)),60,'d',...
               'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],...
               'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);   
           case 3
                   scatter(mean(allMonkeyComparisons(:,1)),mean(allMonkeyComparisons(:,2)),100,'x',...
               'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],...
               'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);   
       end
%        plot(mean(allMonkeyComparisons(:,1)),mean(allMonkeyComparisons(:,2)),'.','MarkerSize',30,...
%            'Color',mcmap(monkeyInd,:));
       monkeyInd = monkeyInd + 1; 
    end
    ax = gca;
    xlim = ax.XLim; ylim = ax.YLim;
    %Diagonal
    maxVal = max(horzcat(xlim,ylim));
    plot([0.5,maxVal],[0.5,maxVal],'--','LineWidth',2,'Color','k')
    %Line at 1
    ax = gca;
    xlim = ax.XLim; ylim = ax.YLim;
    plot([xlim(1),xlim(2)],[1,1],'--','LineWidth',2,'Color','k')
    %Axis Labels
    ylabel('Normalized Error After Shift')
    xlabel('Normalized Error Before Shift')
    
    
    
%% Get overall pct reduction for reporting in text
 monkeyInd = 1;
 for monkey = monkeyList
       allMonkeyComparisons = allMonkeyComparisonsStruct(monkeyInd).allMonkeyComparisons;
       before = mean(allMonkeyComparisons(:,1));
       after = mean(allMonkeyComparisons(:,2));
       monkey
       pctReduction = (((before-1)-(after-1))/(before-1))*100
       monkeyInd = monkeyInd + 1; 
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
    

    
  
    
