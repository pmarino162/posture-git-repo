clear; clc; clf; close all

% Relevant code for computing overt task rxn times:
% Matlab repo / marino / preprocessing / preprocessAndSaveData20220419
% Matlab repo / marino / preprocessing / getKinData20220419

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\Reviewer responses\Analyses\Fig 6';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Datasets to include in analysis 
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};         
    isoDatasetList = {'E20200116','E20200117','E20200120'};    

%% Compute rxn times for all sessions
    resultStruct = struct('monkey',[],'dataset',[],'allRxnTimes',[]);
    structInd = 1;
    task = "bci";
    for datasetList = bciDatasetList
        allRxnTimes = [];
        % Load data
        dataset = datasetList{1,1};
        monkey = dataset(1);
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        % Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',false);  
        [postureList,~,targetList,~,~,~] = getTrajStructDimensions(trajStruct);
        % Get rxnTimes
        if task == "iso" || task == "reach"
            for trial = 1:length(Data)
                rxnTime = Data(trial).kinData.rxnTime;
                allRxnTimes = [allRxnTimes, rxnTime];
            end
        elseif task == "bci"
            for posture = postureList
                for target = targetList
                    allZSmoothFR = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allZSmoothFR;
                    for trial = 1:length(allZSmoothFR)
                        rxnTime = getRxnTimeBCI(target, posture, trial, trajStruct);
                        allRxnTimes = [allRxnTimes,rxnTime];
                    end
                end
            end
        end
        % Store in result struct
        resultStruct(structInd).monkey = monkey;
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).allRxnTimes = allRxnTimes;
        structInd = structInd + 1;
    end
          
%% Combine across monkeys and save
monkeyResultStruct = struct('monkey','','mean',[],'stdErr',[]);
structInd = 1;

for monkey = {'E','N','R'}
    
    thisMonkey = monkey{1};   % extract the character from the cell
    
    % Use strcmp to find entries in resultStruct that match thisMonkey
    monkeyData = resultStruct(strcmp({resultStruct.monkey}, thisMonkey));
    
    
  %  monkeyData = resultStruct([resultStruct.monkey]==monkey);
    combinedRxnTimes = [];
    for i = 1:length(monkeyData)
        combinedRxnTimes = [combinedRxnTimes, monkeyData(i).allRxnTimes];
    end
    
    avgVal = mean(combinedRxnTimes);
    stdVal = std(combinedRxnTimes);
    numObs = length(combinedRxnTimes);
    
    monkeyResultStruct(structInd).monkey = monkey;
    monkeyResultStruct(structInd).mean = avgVal;
    monkeyResultStruct(structInd).stdErr = stdVal / sqrt(numObs);
    
    structInd = structInd + 1;
end


%% Function for computing BCI rxn times

%Walk forwards from start until activity is 50% of the way to
%its peak

function [rxnTime] = getRxnTimeBCI(target, posture, trial, trajStruct)

    % Get trial traj
    allZSmoothFR = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allZSmoothFR;
    time = allZSmoothFR(trial).timestamps;
    traj = allZSmoothFR(trial).traj;
    
    neuralDist = vecnorm(traj,2,2);
    neuralDist = neuralDist - neuralDist(1);
    [maxNeuralDist, maxNeuralDistInd] = max(neuralDist);
    
    threshNeuralDist = 0.5 * maxNeuralDist;

    i = 1;
    while neuralDist(i) < threshNeuralDist
        i = i+1;
        if i == maxNeuralDistInd
           break 
        end
    end
    rxnTime = time(i);
    
end