clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\Reviewer responses\Analyses\Fig 6';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Set parameters
   numPCsToKeep = 10;    %Num PCs to project data into before analysis
   numIterations = 10;
   task = 'Reach'; %Task to analyze ('BCI', 'Iso', or 'Reach')
   statsBinWidth = 0.2; %Difference before shift bin width for stats

   %% Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(5);
    
%% Set up monkey task trial counts struct (determines how many trials to use for each monkey and task)
    trialsToUseStruct = struct('monkey',[],'task',[],'trialsToUse',[]);
    trialsToUseStruct(1) = struct('monkey','E','task','BCI','trialsToUse',12); %These values were computed by taking the minimum # trials in any condition across sessions for each monkey and task (after excluding conditions with less than 10 trials)
    trialsToUseStruct(2) = struct('monkey','N','task','BCI','trialsToUse',28);
    trialsToUseStruct(3) = struct('monkey','R','task','BCI','trialsToUse',16);
    trialsToUseStruct(4) = struct('monkey','E','task','Iso','trialsToUse',38);
    trialsToUseStruct(5) = struct('monkey','E','task','Reach','trialsToUse',14);
    trialsToUseStruct(6) = struct('monkey','N','task','Reach','trialsToUse',32);
    trialsToUseStruct(7) = struct('monkey','R','task','Reach','trialsToUse',44);
    
%% Datasets to include in analysis 
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};         
    isoDatasetList = {'E20200116','E20200117','E20200120'};    

%% Perform all comparisons for all sessions. Save to resultStruct
    if strcmp(task,'BCI')
        taskDatasetList = bciDatasetList;
    elseif strcmp(task,'Iso')
        taskDatasetList = isoDatasetList;
    elseif strcmp(task,'Reach')
        taskDatasetList = reachDatasetList;
    end

    resultStruct = struct('monkey',[],'dataset',[],'result',[]);
    structInd = 1;
    
    for datasetList = taskDatasetList
        %% Set up trajStruct (per session)
        %Load data
        dataset = datasetList{1,1};
        monkey = dataset(1);
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',false);  
        %Project all data down to top PCs
        [allTraj] = collectAllAvgTraj(trajStruct);
        [trajStruct] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep);
        %Remove any conditions for which there weren't enough trials
        trialsToUse = trialsToUseStruct(strcmp({trialsToUseStruct.monkey},monkey) & strcmp({trialsToUseStruct.task},task)).trialsToUse;
        [numCondTrials] = getNumCondTrials(trajStruct,'showHist',false);         
        trajStruct = trajStruct(numCondTrials >= trialsToUse);    
        %get trajStructDims
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
        %Get minimum number of condition trials and timestamps
        [minNumCondTrials] = getMinNumCondTrials(trajStruct);
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
                
        %% Preallocate sessionResultStruct
        sessionResultStruct = struct('predictedPosture',[],'predictedTarget',[],'predictorPosture',[],'predictorTarget',[],... 
            'errorBeforeShift',[],'errorAfterShift',[],'normErrorBeforeShift',[],'normErrorAfterShift',[],...
            'traj1ErrorBeforeShift',[],'traj2ErrorBeforeShift',[],'traj1ErrorAfterShift',[],'traj2ErrorAfterShift',[],'pctChange',[]);
        sessionResultStructInd = 1;
        for predictedPosture = postureList
            withinPostureTargetList = unique([trajStruct([trajStruct.posture]==predictedPosture).target]);       
            for predictedTarget = withinPostureTargetList
                % Get other postures for same target
                for predictorPosture = postureList(postureList>predictedPosture) %In this for loop, predictorTarget=predictedTarget
                    if any([trajStruct.posture]==predictorPosture & [trajStruct.target]==predictedTarget)
                        sessionResultStruct = preallocateSessionResultStructRow(sessionResultStruct, sessionResultStructInd, predictedPosture, predictedTarget, predictorPosture, predictedTarget, numIterations);
                        sessionResultStructInd = sessionResultStructInd + 1;
                    end
                end
                % Get other targets for same posture
                for predictorTarget = withinPostureTargetList(withinPostureTargetList>predictedTarget) %In this loop, predictorPosture=predictedPosture
                    if any([trajStruct.posture]==predictedPosture & [trajStruct.target]==predictorTarget)
                        sessionResultStruct = preallocateSessionResultStructRow(sessionResultStruct, sessionResultStructInd, predictedPosture, predictedTarget, predictedPosture, predictorTarget, numIterations);
                        sessionResultStructInd = sessionResultStructInd + 1;
                    end
                end
            end
        end

        %% Do computation (per session)
        for iteration = 1:numIterations
            %Split trajStruct into two groups
            [trajStruct1,trajStruct2] = splitDataforCV(trajStruct,trialsToUse,binWidth); 
            %Make every comparison
            for i = 1:size(sessionResultStruct,2)
                predictedPosture = sessionResultStruct(i).predictedPosture;
                predictedTarget = sessionResultStruct(i).predictedTarget;
                predictorPosture = sessionResultStruct(i).predictorPosture;
                predictorTarget = sessionResultStruct(i).predictorTarget;
                sessionResultStruct = makeComparisonAndAddToSessionResultStruct(sessionResultStruct, trajStruct1, trajStruct2, predictedPosture, predictedTarget, predictorPosture, predictorTarget, minNumTimestamps, iteration);
            end
        end

        %Compute overall means and pct changes
        for sessionResultStructInd = 1:numel(sessionResultStruct)
           %error
           sessionResultStruct(sessionResultStructInd).meanErrorBeforeShift = mean([sessionResultStruct(sessionResultStructInd).errorBeforeShift]);
           sessionResultStruct(sessionResultStructInd).meanErrorAfterShift = mean([sessionResultStruct(sessionResultStructInd).errorAfterShift]);
           %normalized error 
           sessionResultStruct(sessionResultStructInd).meanNormErrorBeforeShift = mean([sessionResultStruct(sessionResultStructInd).normErrorBeforeShift]);
           sessionResultStruct(sessionResultStructInd).meanNormErrorAfterShift = mean([sessionResultStruct(sessionResultStructInd).normErrorAfterShift]);    
           %pct change
           traj1ErrorBeforeShift =  sessionResultStruct(sessionResultStructInd).traj1ErrorBeforeShift;
           traj2ErrorBeforeShift = sessionResultStruct(sessionResultStructInd).traj2ErrorBeforeShift;
           traj1ErrorAfterShift = sessionResultStruct(sessionResultStructInd).traj1ErrorAfterShift;
           traj2ErrorAfterShift = sessionResultStruct(sessionResultStructInd).traj2ErrorAfterShift;
           errorBeforeShift = sessionResultStruct(sessionResultStructInd).meanErrorBeforeShift;
           errorAfterShift = sessionResultStruct(sessionResultStructInd).meanErrorAfterShift;
           errorBeforeShiftBias = mean(horzcat(traj1ErrorBeforeShift,traj2ErrorBeforeShift));
           errorAfterShiftBias = mean(horzcat(traj1ErrorAfterShift,traj2ErrorAfterShift));
           errorBeforeShiftNoBias = errorBeforeShift - errorBeforeShiftBias;
           errorAfterShiftNoBias = errorAfterShift - errorAfterShiftBias;            
           sessionResultStruct(sessionResultStructInd).pctChange = 100.*(errorAfterShiftNoBias - errorBeforeShiftNoBias)./(errorBeforeShiftNoBias);
        end
        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).monkey = dataset(1);
        resultStruct(structInd).result = sessionResultStruct;
        structInd = structInd + 1;
    end

%% Get Monkey List from resultStruct
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).monkey;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
%% Combine all results across sessions within each monkey. Compute 'diffBetweenGroups' and run t-test
    allMonkeyComparisonStruct = struct('monkey',[],'allMonkeyMeanNormErrorBeforeShift',[],'allMonkeyMeanNormErrorAfterShift',[],'allMonkeyGroup',[],'allMonkeyDiffBetweenGroups',[],'p',[],'h',[],'meanPctChangeAcrossPostures',[],'meanPctChangeAcrossTargets',[]);
    monkeyInd = 1;
    for monkey = monkeyList
        allMonkeyMeanNormErrorBeforeShift = [];
        allMonkeyMeanNormErrorAfterShift = [];
        allMonkeyGroup = [];
        allMonkeyDiffBetweenGroups = [];
        allMonkeyPctChangeAcrossPostures = [];
        allMonkeyPctChangeAcrossTargets = [];
        tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
        for i = 1:size(tempResultStruct,2)
            % Get session result
            result = tempResultStruct(i).result; 
            meanNormErrorBeforeShift = [result.meanNormErrorBeforeShift];
            meanNormErrorAfterShift = [result.meanNormErrorAfterShift];
            pctChange = [result.pctChange];
            
            % Get group numbers 1 = across posture, 2 = smallest target comparison, 0 = all else
            group = zeros(1, length(result));
            group([result.predictedTarget]==[result.predictorTarget]) = 1;
            if ismember(resultStruct(i).dataset, isoDatasetList)
                group(getIndicesOfWithinPostureAngleComparisons(result, 180)) = 2;
            else
                group(getIndicesOfWithinPostureAngleComparisons(result, 45)) = 2;
            end

            % Get before and after shift errors, pct change for each group
            acrossPostureBeforeShift = meanNormErrorBeforeShift(group==1);
            acrossPostureAfterShift = meanNormErrorAfterShift(group==1);
            acrossPosturePctChange = pctChange(group==1);
            acrossTargetBeforeShift = meanNormErrorBeforeShift(group==2);
            acrossTargetAfterShift = meanNormErrorAfterShift(group==2);    
            acrossTargetPctChange = pctChange(group==2);
            
            %Set up bins
            maxBeforeShift = max([acrossPostureBeforeShift, acrossTargetBeforeShift]);
            numBins = ceil(maxBeforeShift/statsBinWidth) + 5;
            binEdges = zeros(1, numBins+1);
            for bin = 1:numBins
                binEdges(bin+1) = binEdges(bin) + statsBinWidth;
            end

            %Bin after shift points based on before shift, take means, get difference between groups, save
            diffBetweenGroups = [];
            for bin = 1:numBins
                binAcrossPostureAfterShift = acrossPostureAfterShift(acrossPostureBeforeShift >= binEdges(bin) & acrossPostureBeforeShift < binEdges(bin + 1));
                binAcrossTargetAfterShift = acrossTargetAfterShift(acrossTargetBeforeShift >= binEdges(bin) & acrossTargetBeforeShift < binEdges(bin + 1));
                if ~isempty(binAcrossPostureAfterShift) && ~isempty(binAcrossTargetAfterShift)
                    diffBetweenGroups = [diffBetweenGroups, (mean(binAcrossTargetAfterShift) - mean(binAcrossPostureAfterShift))];
                end
            end

            % Add all to vectors to combine
            allMonkeyMeanNormErrorBeforeShift = [allMonkeyMeanNormErrorBeforeShift, meanNormErrorBeforeShift];
            allMonkeyMeanNormErrorAfterShift = [allMonkeyMeanNormErrorAfterShift, meanNormErrorAfterShift];
            allMonkeyGroup = [allMonkeyGroup, group];
            allMonkeyDiffBetweenGroups = [allMonkeyDiffBetweenGroups, diffBetweenGroups];
            allMonkeyPctChangeAcrossPostures = [allMonkeyPctChangeAcrossPostures, acrossPosturePctChange];
            allMonkeyPctChangeAcrossTargets = [allMonkeyPctChangeAcrossTargets, acrossTargetPctChange];
        end

    
        % Run stats
        [h, p] = ttest(allMonkeyDiffBetweenGroups);
        
        % Save to struct
        allMonkeyComparisonStruct(monkeyInd).monkey = monkey{1,1};
        allMonkeyComparisonStruct(monkeyInd).allMonkeyMeanNormErrorBeforeShift = allMonkeyMeanNormErrorBeforeShift;
        allMonkeyComparisonStruct(monkeyInd).allMonkeyMeanNormErrorAfterShift = allMonkeyMeanNormErrorAfterShift;
        allMonkeyComparisonStruct(monkeyInd).allMonkeyGroup = allMonkeyGroup;
        allMonkeyComparisonStruct(monkeyInd).allMonkeyDiffBetweenGroups = allMonkeyDiffBetweenGroups;
        allMonkeyComparisonStruct(monkeyInd).h = h;
        allMonkeyComparisonStruct(monkeyInd).p = p;
        allMonkeyComparisonStruct(monkeyInd).meanPctChangeAcrossPostures = mean(allMonkeyPctChangeAcrossPostures);
        allMonkeyComparisonStruct(monkeyInd).meanPctChangeAcrossTargets = mean(allMonkeyPctChangeAcrossTargets);
        monkeyInd = monkeyInd + 1;
    end
    %Save
    %save(fullfile(saveDir,'stats','allMonkeyComparisonStruct.mat'),'allMonkeyComparisonStruct')
    
%% Plot results (per monkey)
for monkey = monkeyList
    figure; hold on;
    % Plot all points
    plotAllPoints(monkey, allMonkeyComparisonStruct, pcmap(1,:), [0 0 0])
    % Plot overall means
    addOverallMeansToPlot(monkey, allMonkeyComparisonStruct, pcmap(1,:), [0 0 0])
    % Format and save
    formatFig6Plot(gca)
    if saveFig
        saveas(gcf,fullfile(saveDir,'scatters',['monkey',monkey{1,1},'_',task,'_scatter.svg']));
    end
end

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

    %Preallocate result struct row
    function sessionResultStruct = preallocateSessionResultStructRow(sessionResultStruct, sessionResultStructInd, predictedPosture, predictedTarget, predictorPosture, predictorTarget, numIterations)
        sessionResultStruct(sessionResultStructInd).predictedPosture = predictedPosture;
        sessionResultStruct(sessionResultStructInd).predictedTarget = predictedTarget;
        sessionResultStruct(sessionResultStructInd).predictorPosture = predictorPosture;
        sessionResultStruct(sessionResultStructInd).predictorTarget = predictorTarget;
        sessionResultStruct(sessionResultStructInd).errorBeforeShift = NaN(1,numIterations);
        sessionResultStruct(sessionResultStructInd).errorAfterShift = NaN(1,numIterations);
        sessionResultStruct(sessionResultStructInd).baselineError = NaN(1,numIterations);
        sessionResultStruct(sessionResultStructInd).normErrorBeforeShift = NaN(1,numIterations);
        sessionResultStruct(sessionResultStructInd).normErrorAfterShift = NaN(1,numIterations);
        sessionResultStruct(sessionResultStructInd).traj1ErrorBeforeShift = NaN(1,numIterations);
        sessionResultStruct(sessionResultStructInd).traj2ErrorBeforeShift = NaN(1,numIterations);
        sessionResultStruct(sessionResultStructInd).traj1ErrorAfterShift = NaN(1,numIterations);
        sessionResultStruct(sessionResultStructInd).traj2ErrorAfterShift = NaN(1,numIterations);
        sessionResultStruct(sessionResultStructInd).pctChange = NaN(1,numIterations);
    end
    
    %Make comparison and add to session result struct
    function sessionResultStruct = makeComparisonAndAddToSessionResultStruct(sessionResultStruct, trajStruct1, trajStruct2, predictedPosture, predictedTarget, predictorPosture, predictorTarget, numPts, iteration)
        predictedCond = find([trajStruct1.target]==predictedTarget & [trajStruct1.posture]==predictedPosture);
        predictorCond = find([trajStruct1.target]==predictorTarget & [trajStruct1.posture]==predictorPosture);
        
        %Get all resampled trajectories
        traj1group1 = trajStruct1(predictedCond).avgPCA.traj;%Predicted traj
        traj1group2 = trajStruct2(predictedCond).avgPCA.traj;
        traj2group1 = trajStruct1(predictorCond).avgPCA.traj;%Comparison traj
        traj2group2 = trajStruct2(predictorCond).avgPCA.traj;

        %Compute traj1 diff before shift (w/in condition)
        tempNumPts = getTempNumPts(traj1group1,traj1group2,numPts);
        traj1ErrorBeforeShift = getMeanDist(traj1group1,traj1group2,tempNumPts);

        %Compute traj1 diff after shift (w/in %condition)
        alpha = getOptimalAlpha(traj1group1,traj1group2,tempNumPts);
        shiftTraj1Group2 = traj1group2 + alpha;  
        traj1ErrorAfterShift = getMeanDist(traj1group1,shiftTraj1Group2,tempNumPts);

        %Compute traj2 diff before shift (w/in condition)
        tempNumPts = getTempNumPts(traj2group1,traj2group2,numPts);
        traj2ErrorBeforeShift = getMeanDist(traj2group1,traj2group2,tempNumPts);

        %Compute traj2 diff after shift (w/in %condition)
        alpha = getOptimalAlpha(traj2group1,traj2group2,tempNumPts);
        shiftTraj2Group2 = traj2group2 + alpha;  
        traj2ErrorAfterShift = getMeanDist(traj2group1,shiftTraj2Group2,tempNumPts);

        %Compute diff before shift
        tempNumPts = getTempNumPts(traj1group1,traj2group1,numPts);
        errorBeforeShift = getMeanDist(traj1group1,traj2group1,tempNumPts);
        normErrorBeforeShift = errorBeforeShift./traj1ErrorAfterShift;

        %Compute diff after shift
        alpha = getOptimalAlpha(traj1group1,traj2group1,tempNumPts);
        shiftTraj2Group1 = traj2group1 + alpha;  
        errorAfterShift = getMeanDist(traj1group1,shiftTraj2Group1,tempNumPts);
        normErrorAfterShift = errorAfterShift./traj1ErrorAfterShift;

        %Store results in sessionResultStruct
        sessionResultStructInd = find([sessionResultStruct.predictedPosture]==predictedPosture & ...
            [sessionResultStruct.predictedTarget]==predictedTarget & [sessionResultStruct.predictorPosture]==predictorPosture & ...
            [sessionResultStruct.predictorTarget]==predictorTarget);
        sessionResultStruct(sessionResultStructInd).errorBeforeShift(iteration) = errorBeforeShift;
        sessionResultStruct(sessionResultStructInd).errorAfterShift(iteration) = errorAfterShift;     
        sessionResultStruct(sessionResultStructInd).normErrorBeforeShift(iteration) = normErrorBeforeShift;
        sessionResultStruct(sessionResultStructInd).normErrorAfterShift(iteration) = normErrorAfterShift;
        sessionResultStruct(sessionResultStructInd).traj1ErrorBeforeShift(iteration) = traj1ErrorBeforeShift;
        sessionResultStruct(sessionResultStructInd).traj2ErrorBeforeShift(iteration) = traj2ErrorBeforeShift;
        sessionResultStruct(sessionResultStructInd).traj1ErrorAfterShift(iteration) = traj1ErrorAfterShift;
        sessionResultStruct(sessionResultStructInd).traj2ErrorAfterShift(iteration) = traj2ErrorAfterShift;
    end
    
    %Get indices of within posture target comparisons of particular angle
    function indices = getIndicesOfWithinPostureAngleComparisons(result, angle)
        % Unpack
        predictedPosture = [result.predictedPosture];
        predictedTarget = [result.predictedTarget];
        predictorPosture = [result.predictorPosture];
        predictorTarget = [result.predictorTarget];

        % Convert to angle
        predictedTargetAngles = 45*(predictedTarget-1);
        predictorTargetAngles = 45*(predictorTarget-1);
        
        % Compute the circular difference in target angles
        diffAngles = mod(abs(predictedTargetAngles - predictorTargetAngles), 360);

        % Adjust differences greater than 180
        diffAngles(diffAngles > 180) = 360 - diffAngles(diffAngles > 180);

        % Find indices where postures are equal and angle difference matches the input
        indices = find(predictedPosture == predictorPosture & diffAngles == angle);  
    end
    
    %Get number of points used in comparison
    function tempNumPts = getTempNumPts(traj1,traj2,numPts)
        if min([size(traj1,1),size(traj2,1)]) < numPts
            tempNumPts = min([size(traj1,1),size(traj2,1)]);
        else
            tempNumPts = numPts;
        end
    end
    
    %Plot all points
    function plotAllPoints(monkey, allMonkeyComparisonStruct, acrossPostureColor, acrossTargetColor)
        tempComparisonStruct = allMonkeyComparisonStruct(strcmp({allMonkeyComparisonStruct.monkey},monkey{1,1}));
        meanNormErrorBeforeShift = tempComparisonStruct.allMonkeyMeanNormErrorBeforeShift;
        meanNormErrorAfterShift = tempComparisonStruct.allMonkeyMeanNormErrorAfterShift;
        group = tempComparisonStruct.allMonkeyGroup;
        
       scatter(meanNormErrorBeforeShift(group==1),meanNormErrorAfterShift(group==1),25,'o',...
           'MarkerEdgeColor', acrossPostureColor, 'MarkerFaceColor', acrossPostureColor,...
           'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',.5);   

       scatter(meanNormErrorBeforeShift(group==2),meanNormErrorAfterShift(group==2),25,'o',...
           'MarkerEdgeColor', acrossTargetColor, 'MarkerFaceColor', acrossTargetColor,...
           'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',.5); 
    end
    
    %Add overall means
    function addOverallMeansToPlot(monkey, allMonkeyComparisonStruct, acrossPostureColor, acrossTargetColor)
        tempComparisonStruct = allMonkeyComparisonStruct(strcmp({allMonkeyComparisonStruct.monkey},monkey{1,1}));
        meanNormErrorBeforeShift = tempComparisonStruct.allMonkeyMeanNormErrorBeforeShift;
        meanNormErrorAfterShift = tempComparisonStruct.allMonkeyMeanNormErrorAfterShift;
        group = tempComparisonStruct.allMonkeyGroup;
    
       scatter(mean(meanNormErrorBeforeShift(group==1)),mean(meanNormErrorAfterShift(group==1)),275,'o',...
            'MarkerEdgeColor', acrossPostureColor, 'MarkerFaceColor',[1 1 1],...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);  

        scatter(mean(meanNormErrorBeforeShift(group==2)),mean(meanNormErrorAfterShift(group==2)),275,'o',...
            'MarkerEdgeColor', acrossTargetColor, 'MarkerFaceColor',[1 1 1],...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);  
    end
    
    %Format plot
    function formatFig6Plot(ax)
        
        xlim = ax.XLim; ylim = ax.YLim;
        %Diagonal Line
        maxVal = 7;
        minVal = 0.75;
        plot([minVal ,maxVal],[minVal, maxVal],'--','LineWidth',2,'Color','k')
        %Line at 1
        xlim = ax.XLim; ylim = ax.YLim;
        plot([xlim(1),xlim(2)],[1,1],'--','LineWidth',2,'Color','k')
        %Axis Labels, ticks
        ax.XLim = [minVal, maxVal];
        ax.YLim = [minVal, maxVal];
        xticks([1,xlim(2)])
        yticks([1,ylim(2)])
        ax.TickDir = 'out';
        set(ax,'fontname','arial'); set(gca,'fontsize',20);
    end