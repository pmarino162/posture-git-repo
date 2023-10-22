clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 6';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set parameters
   numPCsToKeep = 10;    %Num PCs to project data into before analysis
   numIterations = 1000;
   cutoffNumTraj = 10; %Num trials that must be present in a condition to keep it for analysis 
   
%% Datasets to include in analysis 
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};         
    isoDatasetList = {'E20200116','E20200117','E20200120'};    

%% Main loop
    resultStruct = struct('animal',[],'dataset',[],'result',[]);
    structInd = 1;
    for datasetList = reachDatasetList%{'E20200316'}%{'N20190222','N20190226','R20200221','R20200222'}%reachDatasetList%{'E20200316'}%bciDatasetList% reachDatasetList%{'E20210707','N20190226','R20200221'}%{'E20200316','N20171215','R20201020'}%{ 'R20200221'}%bciDatasetList%{'E20210706','E20210707','E20210708','E20210709'}%reachDatasetList
        %% Set up trajStruct
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',false);  
        %Remove any conditions for which there weren't enough trials
        [numCondTrials] = getNumCondTrials(trajStruct,'showHist',false);         
        trajStruct = trajStruct(numCondTrials >= cutoffNumTraj);    
        %get trajStructDims
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
        %Get minimum number of condition trials and timestamps
        [minNumCondTrials] = getMinNumCondTrials(trajStruct);
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        %Project all data down to top PCs
        [allTraj] = collectAllAvgTraj(trajStruct);
        [trajStruct] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep);
                
        %% Preallocate sessionResultStruct
        sessionResultStruct = struct('predictedPosture',[],'predictedTarget',[],'predictorPosture',[],'predictorTarget',[],... 
            'errorBeforeShift',[],'errorAfterShift',[],'normErrorBeforeShift',[],'normErrorAfterShift',[],...
            'R2BeforeShift',[],'R2AfterShift',[],'pctChange',[]);
        sessionResultStructInd = 1;
        for predictedPosture = postureList
            withinPostureTargetList = unique([trajStruct([trajStruct.posture]==predictedPosture).target]);       
            for predictedTarget = withinPostureTargetList
                for predictorPosture = postureList(postureList>predictedPosture)
                    if any([trajStruct.posture]==predictorPosture & [trajStruct.target]==predictedTarget)
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
                        sessionResultStruct(sessionResultStructInd).pctChange = NaN(1,numIterations);
                        sessionResultStructInd = sessionResultStructInd + 1;
                    end
                end
            end
        end
        
        %% Do computation
        numPts = minNumTimestamps;
        for i = 1:numIterations
            i
            %Split trajStruct into two groups
            [trajStruct1,trajStruct2] = splitDataforCV(trajStruct,minNumCondTrials,binWidth); 
            %Make every comparison
            for predictedPosture = postureList
                withinPostureTargetList = unique([trajStruct([trajStruct.posture]==predictedPosture).target]);       
                for predictedTarget = withinPostureTargetList
                    for predictorPosture = postureList(postureList>predictedPosture)
                        if any([trajStruct.posture]==predictorPosture & [trajStruct.target]==predictedTarget)
                            predictedCond = find([trajStruct1.target]==predictedTarget & [trajStruct1.posture]==predictedPosture);
                            predictorCond = find([trajStruct2.target]==predictedTarget & [trajStruct2.posture]==predictorPosture);
                            %Compute comparison difference 
                            traj1 = trajStruct1(predictedCond).avgPCA.traj;%Predicted traj
                            traj2 = trajStruct2(predictorCond).avgPCA.traj;%Comparison traj
                            baselineTraj = trajStruct2(predictedCond).avgPCA.traj;%Used for normalization
                            unshiftTraj2 = traj2;
                            if min([size(traj1,1),size(traj2,1)]) < numPts
                                tempNumPts = min([size(traj1,1),size(traj2,1)]);
                            else
                                tempNumPts = numPts;
                            end
                            alpha = getOptimalAlpha(traj1,traj2,tempNumPts);
                            shiftTraj2 = traj2 + alpha;  
                            alpha = getOptimalAlpha(traj1,baselineTraj,tempNumPts);
                            shiftBaselineTraj = baselineTraj + alpha;
                            %Compute absolute and normalized errors
                            errorBeforeShift = getMeanDist(traj1,unshiftTraj2,tempNumPts);
                            errorAfterShift = getMeanDist(traj1,shiftTraj2,tempNumPts);
                            baselineError = getMeanDist(traj1,shiftBaselineTraj,tempNumPts);
                            normErrorBeforeShift = errorBeforeShift./baselineError;
                            normErrorAfterShift = errorAfterShift./baselineError;
                            R2BeforeShift = getR2(traj1,unshiftTraj2,tempNumPts);
                            R2AfterShift = getR2(traj1,shiftTraj2,tempNumPts);
                            %Store results in sessionResultStruct
                            sessionResultStructInd = find([sessionResultStruct.predictedPosture]==predictedPosture & ...
                                [sessionResultStruct.predictedTarget]==predictedTarget & [sessionResultStruct.predictorPosture]==predictorPosture);
                            sessionResultStruct(sessionResultStructInd).errorBeforeShift(i) = errorBeforeShift;
                            sessionResultStruct(sessionResultStructInd).errorAfterShift(i) = errorAfterShift;     
                            sessionResultStruct(sessionResultStructInd).normErrorBeforeShift(i) = normErrorBeforeShift;
                            sessionResultStruct(sessionResultStructInd).normErrorAfterShift(i) = normErrorAfterShift;
                            sessionResultStruct(sessionResultStructInd).R2BeforeShift(i) = R2BeforeShift;
                            sessionResultStruct(sessionResultStructInd).R2AfterShift(i) = R2AfterShift;
                            sessionResultStruct(sessionResultStructInd).pctChange(i) = 100*(((errorAfterShift-baselineError)-(errorBeforeShift-baselineError))/(errorBeforeShift-baselineError));
                        end
                    end
                end
            end
        end
           
        %Compute overall means
        for sessionResultStructInd = 1:numel(sessionResultStruct)
           %R2
           sessionResultStruct(sessionResultStructInd).meanR2BeforeShift = mean([sessionResultStruct(sessionResultStructInd).R2BeforeShift]);
           sessionResultStruct(sessionResultStructInd).meanR2AfterShift = mean([sessionResultStruct(sessionResultStructInd).R2AfterShift]);
           %error
           sessionResultStruct(sessionResultStructInd).meanErrorBeforeShift = mean([sessionResultStruct(sessionResultStructInd).errorBeforeShift]);
           sessionResultStruct(sessionResultStructInd).meanErrorAfterShift = mean([sessionResultStruct(sessionResultStructInd).errorAfterShift]);
           %normalized error 
           sessionResultStruct(sessionResultStructInd).meanNormErrorBeforeShift = mean([sessionResultStruct(sessionResultStructInd).normErrorBeforeShift]);
           sessionResultStruct(sessionResultStructInd).meanNormErrorAfterShift = mean([sessionResultStruct(sessionResultStructInd).normErrorAfterShift]);    
           %pctChange
           sessionResultStruct(sessionResultStructInd).meanPctChange = mean([sessionResultStruct(sessionResultStructInd).pctChange]);
        end
       
        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).result = sessionResultStruct;
        structInd = structInd + 1;
    end

%% Plot results
    %Get monkeyList
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
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
               normErrorBeforeShiftMean = result(comparison).meanNormErrorBeforeShift;
               normErrorAfterShiftMean = result(comparison).meanNormErrorAfterShift;
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
    
    %Plot overall means
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
       monkeyInd = monkeyInd + 1; 
    end
    ax = gca;
    xlim = ax.XLim; ylim = ax.YLim;
    
    %Diagonal Line
    maxVal = max(horzcat(xlim,ylim));
    plot([0.5,maxVal],[0.5,maxVal],'--','LineWidth',2,'Color','k')
    %Line at 1
    ax = gca;
    xlim = ax.XLim; ylim = ax.YLim;
    plot([xlim(1),xlim(2)],[1,1],'--','LineWidth',2,'Color','k')
    %Axis Labels, ticks
    xticks([1,xlim(2)])
    yticks([1,ylim(2)])
    ax.TickDir = 'out';
    ylabel('Normalized Error After Shift')
    xlabel('Normalized Error Before Shift')
    
%% Get overall pct reduction for reporting in text
    for monkey = monkeyList
        monkey
        tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
        pctChangeDist = [];
        for session = 1:numel(tempResultStruct)
            result = tempResultStruct(session).result;
            %pctChangeDist = [pctChangeDist,result.pctChange];
            pctChangeDist = [pctChangeDist,result.meanPctChange];
        end
        meanPctReduction = mean(pctChangeDist)
        SEMPctReduction = std(pctChangeDist)/sqrt(length(pctChangeDist))
        figure
        histogram(pctChangeDist)
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
    

    
  
    
