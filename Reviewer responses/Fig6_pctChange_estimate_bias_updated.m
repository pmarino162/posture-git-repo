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
        dataset
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',false);  
        %Project all data down to top PCs
        [allTraj] = collectAllAvgTraj(trajStruct);
        [trajStruct] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep);
        %Remove any conditions for which there weren't enough trials
        [numCondTrials] = getNumCondTrials(trajStruct,'showHist',false);         
        trajStruct = trajStruct(numCondTrials >= cutoffNumTraj);    
        %get trajStructDims
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
        %Get minimum number of condition trials and timestamps
        [minNumCondTrials] = getMinNumCondTrials(trajStruct);
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        
        %% Preallocate sessionResultStruct
        sessionResultStruct = struct('predictedPosture',[],'predictedTarget',[],'predictorPosture',[],'predictorTarget',[],... 
            'diffBeforeShift',[],'diffAfterShift',[],'traj1DiffBeforeShift',[],'traj2DiffBeforeShift',[],...
            'traj1DiffAfterShift',[],'traj2DiffAfterShift',[],'pctChange',[]);
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
                        sessionResultStruct(sessionResultStructInd).diffBeforeShift = NaN(1,numIterations);
                        sessionResultStruct(sessionResultStructInd).diffAfterShift = NaN(1,numIterations);
                        sessionResultStruct(sessionResultStructInd).traj1DiffBeforeShift = NaN(1,numIterations);
                        sessionResultStruct(sessionResultStructInd).traj2DiffBeforeShift = NaN(1,numIterations);
                        sessionResultStruct(sessionResultStructInd).traj1DiffAfterShift = NaN(1,numIterations);
                        sessionResultStruct(sessionResultStructInd).traj2DiffAfterShift = NaN(1,numIterations);
                        sessionResultStruct(sessionResultStructInd).pctChange = NaN;
                        sessionResultStructInd = sessionResultStructInd + 1;
                    end
                end
            end
        end
        
         
        %% Do computation
        numPts = minNumTimestamps;
        for i = 1:numIterations
            %Split trajStruct into two groups
            [trajStruct1,trajStruct2] = splitDataforCV(trajStruct,minNumCondTrials,binWidth); 
            %Make every comparison
            for predictedPosture = postureList
                withinPostureTargetList = unique([trajStruct([trajStruct.posture]==predictedPosture).target]);       
                for predictedTarget = withinPostureTargetList
                    for predictorPosture = postureList(postureList>predictedPosture)
                        if any([trajStruct.posture]==predictorPosture & [trajStruct.target]==predictedTarget)
                            predictedCond = find([trajStruct1.target]==predictedTarget & [trajStruct1.posture]==predictedPosture);
                            predictorCond = find([trajStruct1.target]==predictedTarget & [trajStruct1.posture]==predictorPosture);
                            
                            %Get all resampled avg trajectories
                            traj1group1 = trajStruct1(predictedCond).avgPCA.traj;%Predicted traj
                            traj1group2 = trajStruct2(predictedCond).avgPCA.traj;
                            traj2group1 = trajStruct1(predictorCond).avgPCA.traj;%Comparison traj
                            traj2group2 = trajStruct2(predictorCond).avgPCA.traj;
                            
                            %Compute diff before shift
                            tempNumPts = getTempNumPts(traj1group1,traj2group1,numPts);
                            diffBeforeShift = getMeanDist(traj1group1,traj2group1,tempNumPts);
                            
                            %Compute diff after shift
                            alpha = getOptimalAlpha(traj1group1,traj2group1,tempNumPts);
                            shiftTraj2Group1 = traj2group1 + alpha;  
                            diffAfterShift = getMeanDist(traj1group1,shiftTraj2Group1,tempNumPts);
                            
                            %Compute traj1 diff before shift (w/in condition)
                            tempNumPts = getTempNumPts(traj1group1,traj1group2,numPts);
                            traj1DiffBeforeShift = getMeanDist(traj1group1,traj1group2,tempNumPts);
                            
                            %Compute traj1 diff after shift (w/in %condition)
                            alpha = getOptimalAlpha(traj1group1,traj1group2,tempNumPts);
                            shiftTraj1Group2 = traj1group2 + alpha;  
                            traj1DiffAfterShift = getMeanDist(traj1group1,shiftTraj1Group2,tempNumPts);
                            
                            %Compute traj2 diff before shift (w/in condition)
                            tempNumPts = getTempNumPts(traj2group1,traj2group2,numPts);
                            traj2DiffBeforeShift = getMeanDist(traj2group1,traj2group2,tempNumPts);
                            
                            %Compute traj2 diff after shift (w/in %condition)
                            alpha = getOptimalAlpha(traj2group1,traj2group2,tempNumPts);
                            shiftTraj2Group2 = traj2group2 + alpha;  
                            traj2DiffAfterShift = getMeanDist(traj2group1,shiftTraj2Group2,tempNumPts);
                            

                            %Store results in sessionResultStruct
                            sessionResultStructInd = find([sessionResultStruct.predictedPosture]==predictedPosture & ...
                                [sessionResultStruct.predictedTarget]==predictedTarget & [sessionResultStruct.predictorPosture]==predictorPosture);
                            sessionResultStruct(sessionResultStructInd).diffBeforeShift(i) = diffBeforeShift;
                            sessionResultStruct(sessionResultStructInd).diffAfterShift(i) = diffAfterShift;     
                            sessionResultStruct(sessionResultStructInd).traj1DiffBeforeShift(i) = traj1DiffBeforeShift;
                            sessionResultStruct(sessionResultStructInd).traj2DiffBeforeShift(i) = traj2DiffBeforeShift;
                            sessionResultStruct(sessionResultStructInd).traj1DiffAfterShift(i) = traj1DiffAfterShift;
                            sessionResultStruct(sessionResultStructInd).traj2DiffAfterShift(i) = traj2DiffAfterShift;
                            
                        end
                    end
                end
            end
        end
        
        %Compute percent changes
        for sessionResultStructInd = 1:numel(sessionResultStruct)                       
            diffBeforeShift = mean([sessionResultStruct(sessionResultStructInd).diffBeforeShift]);
            diffAfterShift = mean([sessionResultStruct(sessionResultStructInd).diffAfterShift]);
            traj1DiffBeforeShift =  sessionResultStruct(sessionResultStructInd).traj1DiffBeforeShift;
            traj2DiffBeforeShift = sessionResultStruct(sessionResultStructInd).traj2DiffBeforeShift;
            traj1DiffAfterShift = sessionResultStruct(sessionResultStructInd).traj1DiffAfterShift;
            traj2DiffAfterShift = sessionResultStruct(sessionResultStructInd).traj2DiffAfterShift;
            diffBeforeShiftBias = mean(horzcat(traj1DiffBeforeShift,traj2DiffBeforeShift));
            diffAfterShiftBias = mean(horzcat(traj1DiffAfterShift,traj2DiffAfterShift));
            diffBeforeShiftNoBias = diffBeforeShift - diffBeforeShiftBias;
            diffAfterShiftNoBias = diffAfterShift - diffAfterShiftBias;
            sessionResultStruct(sessionResultStructInd).pctChange = 100.*(diffAfterShiftNoBias - diffBeforeShiftNoBias)./(diffBeforeShiftNoBias);
        end
        
        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).result = sessionResultStruct;
        structInd = structInd + 1;
    end
    
%% Get overall pct reduction for reporting in text
    %Get monkeyList
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
    for monkey = monkeyList
        monkey
        tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
        pctChangeDist = [];
        for session = 1:numel(tempResultStruct)
            result = tempResultStruct(session).result;
            pctChangeDist = [pctChangeDist,result.pctChange];
        end
        meanPctReduction = mean(pctChangeDist)
        SEMPctReduction = std(pctChangeDist)/sqrt(length(pctChangeDist))
        figure
        histogram(pctChangeDist)
    end    
    
    
%% Local functions
    function tempNumPts = getTempNumPts(traj1,traj2,numPts)
        if min([size(traj1,1),size(traj2,1)]) < numPts
            tempNumPts = min([size(traj1,1),size(traj2,1)]);
        else
            tempNumPts = numPts;
        end
    end

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