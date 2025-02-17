clear; clc;

%% Load data from session
    %EARL Data = {'E20200316','E20200317','E20200318','E20200319'};
    dataset = 'E20200318';   
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    restTrialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
    restTrialInclStates(1).inclStates = {{'state','Touch Bar Reach','first',50},{'state','Step 1','first',0}};
    
    %NIGEL/ DATA {'N20171215','N20180221','R20201020','R20201021'}; 
%     dataset = 'N20171215';   
%     
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     
%     trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
%     trialInclStates(1).inclStates = {{'state','Cursor Release','first',0},{'state','Target Hold','first',0}};
%     
%     restTrialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
%     restTrialInclStates(1).inclStates = {{'state','Start','first',25},{'state','Center Exit','first',0}};
    
%     %ROCKY
%     dataset = 'R20201020';   
%     
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     
%     trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
%     trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
%     
%     restTrialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
%     restTrialInclStates(1).inclStates = {{'state','Center','first',25},{'state','React','first',0}};
    
    
    % Common to all
    [Data,zScoreParams] = loadData(dataset);
    [Data] = removeShortBCIandIsoTrials(Data,dataset);
    binWidth = 25;kernelStdDev = 25; trajFields = {'zSmoothFR'};
    
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams); 
    [postureList,~,targetList,~,~,~] = getTrajStructDimensions(trajStruct);
    baselineTrajStruct = getTrajStruct(Data,condFields,trajFields,restTrialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);
    
    
    
%% Plot neural distance for a trial
    target = 3;
    posture = 1;
    trial = 5;
    
    rxnTime = getRxnTime(target, posture, trial, trajStruct, baselineTrajStruct, true);

%% Compute for all trials and produce distribution
    allRxnTimes = [];
    for posture = postureList
        for target = targetList
            allZSmoothFR = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allZSmoothFR;
            for trial = 1:length(allZSmoothFR)
                rxnTime = getRxnTime(target, posture, trial, trajStruct, baselineTrajStruct, false);
                allRxnTimes = [allRxnTimes,rxnTime];
            end
        end
    end
    figure;
    histogram(allRxnTimes);
    xlabel('rxn time (ms)');
    ylabel('count');

%% Function for computing reaction times

function [rxnTime] = getRxnTime(target, posture, trial, trajStruct, baselineTrajStruct, plotTrial)

    allZSmoothFR = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allZSmoothFR;
    allZSmoothFRBaseline = baselineTrajStruct([baselineTrajStruct.posture]==posture & [baselineTrajStruct.target]==target).allZSmoothFR;
    
    time = allZSmoothFR(trial).timestamps;
    trialNum = allZSmoothFR(trial).trialNum;
    traj = allZSmoothFR(trial).traj;
    
    timeBaseline = allZSmoothFRBaseline([allZSmoothFRBaseline.trialNum]==trialNum).timestamps;  
    trajBaseline = allZSmoothFRBaseline([allZSmoothFRBaseline.trialNum]==trialNum).traj;

    
    
    baseline = mean(trajBaseline,1);
    
    %neuralDist = vecnorm(traj-baseline,2,2);
    
    %restNeuralDist = vecnorm(trajBaseline-baseline,2,2);
    
    
    neuralDist = vecnorm(traj,2,2);
    
    restNeuralDist = vecnorm(trajBaseline,2,2);
    
    
    %Algo 1: walk forwards from start until activity is 50% of the way to
    %its peak
    
    neuralDist = neuralDist - neuralDist(1);
    
    [maxNeuralDist, maxNeuralDistInd] = max(neuralDist);
    
    threshNeuralDist = 0.5 * maxNeuralDist;
    % Step backwards to get 20% mark
    i = 1;
    while neuralDist(i) < threshNeuralDist
        i = i+1;
        if i == maxNeuralDistInd
           break 
        end
    end
    rxnTime = time(i);
    
    if plotTrial == true
        figure
        subplot(1,2,2); 
        hold on;
        plot(time, neuralDist);
        ax2 = gca;
        xlimits = ax2.XLim; 
        yLimits2 = ax2.YLim;       % store the y-limits after plotting

        plot(xlimits, threshNeuralDist*ones(1,2),'--r')
        plot(rxnTime*ones(1,2), yLimits2, '-r')

        xlabel('time (ms)');
        ylabel('neural dist');

        % --- Subplot #1 ---
        subplot(1,2,1);
        plot(timeBaseline, restNeuralDist);
        ax1 = gca;
        yLimits1 = ax1.YLim;       % store the y-limits after plotting

        xlabel('time (ms)');
        ylabel('neural dist');

        % --- Compute common y-limits ---
        yMin = min([yLimits1(1), yLimits2(1)]);
        yMax = max([yLimits1(2), yLimits2(2)]);

        % --- Apply new y-limits to both subplots ---
        ax1.YLim = [yMin, yMax];
        ax2.YLim = [yMin, yMax];

    end
    
    

end
