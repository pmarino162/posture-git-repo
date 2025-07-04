clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 6';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set parameters
   numPCsToKeep = 10;    %Num PCs to project data into before analysis
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
    for datasetList = isoDatasetList%{'E20200316'}%{'N20190222','N20190226','R20200221','R20200222'}%reachDatasetList%{'E20200316'}%bciDatasetList% reachDatasetList%{'E20210707','N20190226','R20200221'}%{'E20200316','N20171215','R20201020'}%{ 'R20200221'}%bciDatasetList%{'E20210706','E20210707','E20210708','E20210709'}%reachDatasetList
        %% Set up trajStruct
        %Load data
        dataset = datasetList{1,1};
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

                
        %% Compute sessionResultStruct
        numPts = minNumTimestamps;
        sessionResultStruct = struct('predictedPosture',[],'predictedTarget',[],'predictorPosture',[],'predictorTarget',[],... 
            'errorBeforeShift',[],'errorAfterShift',[],'normErrorBeforeShift',[],'normErrorAfterShift',[],'pctChange',[]);
        sessionResultStructInd = 1;
        for predictedPosture = postureList
            withinPostureTargetList = unique([trajStruct([trajStruct.posture]==predictedPosture).target]);       
            for predictedTarget = withinPostureTargetList
                for predictorPosture = postureList(postureList>predictedPosture)
                    if any([trajStruct.posture]==predictorPosture & [trajStruct.target]==predictedTarget)

                        predictedCond = find([trajStruct.target]==predictedTarget & [trajStruct.posture]==predictedPosture);
                        predictorCond = find([trajStruct.target]==predictedTarget & [trajStruct.posture]==predictorPosture);
                        
                        %Compute comparison difference 
                        traj1 = trajStruct(predictedCond).avgPCA.traj;%Predicted traj
                        traj2 = trajStruct(predictorCond).avgPCA.traj;%Comparison traj

                        unshiftTraj2 = traj2;
                        if min([size(traj1,1),size(traj2,1)]) < numPts
                            tempNumPts = min([size(traj1,1),size(traj2,1)]);
                        else
                            tempNumPts = numPts;
                        end
                        
                        alpha = getOptimalAlpha(traj1,traj2,tempNumPts);
                        shiftTraj2 = traj2 + alpha;  

                        %Compute absolute and normalized errors
                        errorBeforeShift = getMeanDist(traj1,unshiftTraj2,tempNumPts);
                        errorAfterShift = getMeanDist(traj1,shiftTraj2,tempNumPts);

                        sessionResultStruct(sessionResultStructInd).predictedPosture = predictedPosture;
                        sessionResultStruct(sessionResultStructInd).predictedTarget = predictedTarget;
                        sessionResultStruct(sessionResultStructInd).predictorPosture = predictorPosture;
                        sessionResultStruct(sessionResultStructInd).predictorTarget = predictedTarget;
                        sessionResultStruct(sessionResultStructInd).errorBeforeShift = errorBeforeShift;
                        sessionResultStruct(sessionResultStructInd).errorAfterShift = errorAfterShift;
                        sessionResultStruct(sessionResultStructInd).pctChange = 100*(errorAfterShift - errorBeforeShift)/errorBeforeShift;
                        sessionResultStructInd = sessionResultStructInd + 1;
                    end
                end
            end
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
            %pctChangeDist = [pctChangeDist,result.meanPctChange];
            %pc
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
    
    
  
    
