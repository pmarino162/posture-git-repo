clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Figure 2';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Load data
    dataset = 'E20200317';
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
        %Reach
        case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
            trialInclStates(1).trialName = {'GridReaching'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
        case {'R20200221','R20200222'}
            trialInclStates(1).trialName = {'Rocky Dissociation'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
        case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
            trialInclStates(1).trialName = {'Nigel Dissociation'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    %Eliminate postures for which there are not trials to all targets
    switch dataset
        case {'E20210706'}
            trajStruct([trajStruct.posture]==5 | [trajStruct.posture]==6) = [];
    end
    
%% Get timestamps dist and min number timestamps
    numTimestamps = [];
    numCondTraj = [];
    for i = 1:size(trajStruct,2)
       numTraj = size(trajStruct(i).allSmoothFR,2);
       numCondTraj = [numCondTraj,numTraj];
       for j = 1:numTraj
          numTimestamps = [numTimestamps,size(trajStruct(i).allSmoothFR(j).timestamps,2)]; 
       end
    end
    figure
    histogram(numTimestamps)
    xlabel('Number of 25ms bins')
    ylabel('Number of trials')
    
    figure
    histogram(numCondTraj)
    xlabel('Number of trials')
    ylabel('Number of conditions')
    
% Get minimum number of trials and timestamps
    [minNumTimestamps,i] = min(numTimestamps);
    [minNumCondTraj,i] = min(numCondTraj);
    
    %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);

%% Do computation
    %Parameters
    numIterations = 100;
    numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw
    numPts = minNumTimestamps; %Number of points from each trajectory

    %Get lower bound and error reduction 
    pctError = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
    index = 1;
    for target = targetList %for every target, make every posture comparison
       for posture1Ind = 1:numPostures-1
           posture1Ind
           posture1 = postureList(posture1Ind);
           for posture2Ind = posture1Ind+1:numPostures  
               posture2 = postureList(posture2Ind);
               cond1 = find([trajStruct.target]==target & [trajStruct.posture]==posture1);
               numTraj1 = size(trajStruct(cond1).allSmoothFR,2);
               cond2 = find([trajStruct.target]==target & [trajStruct.posture]==posture2);
               numTraj2 = size(trajStruct(cond2).allSmoothFR,2);
               for i = 1:numIterations
                    %Select half the trials from each condition 
                    sampInd1 = randsample(numTraj1,numSample);
                    traj1Struct = trajStruct(cond1).allSmoothFR(sampInd1);
                    sampInd2 = randsample(numTraj2,numSample);
                    traj2Struct = trajStruct(cond2).allSmoothFR(sampInd2);
                    %Form condition averages
                    traj1 = getAvgTraj20211210(traj1Struct,binWidth);
                    traj2 = getAvgTraj20211210(traj2Struct,binWidth);
                    %Get average error between trajectories
                    totalError = getMeanDist(traj1,traj2,numPts);
                    %Subtract mean from group 2
                    alpha = getOptimalAlpha(traj1,traj2,numPts);
                    traj2 = traj2 + alpha;
                    %Get remaining average error
                    shiftError = getMeanDist(traj1,traj2,numPts);
                    %Perform self-comparison
                    selfCompCond = randi(2);
                    if selfCompCond == 1
                        numTrajRemaining = numTraj1 - numSample;    
                        selfSampInd = randsample(numTrajRemaining,numSample);
                        remainingInd = setdiff(1:numTraj1,sampInd1);
                        selfSampInd = remainingInd(selfSampInd);
                        selfTrajStruct = trajStruct(cond1).allSmoothFR(selfSampInd);
                        selfTraj = getAvgTraj20211210(selfTrajStruct,binWidth);
                        alpha = getOptimalAlpha(traj1,selfTraj,numPts);
                        selfTraj = selfTraj + alpha;
                        lowerBound = getMeanDist(traj1,selfTraj,numPts);
                    elseif selfCompCond == 2
                        numTrajRemaining = numTraj2 - numSample;   
                        selfSampInd = randsample(numTrajRemaining,numSample);
                        remainingInd = setdiff(1:numTraj2,sampInd2);
                        selfSampInd = remainingInd(selfSampInd);
                        selfTrajStruct = trajStruct(cond2).allSmoothFR(selfSampInd);
                        selfTraj = getAvgTraj20211210(selfTrajStruct,binWidth);
                        alpha = getOptimalAlpha(traj2,selfTraj,numPts);
                        selfTraj = selfTraj + alpha;
                        lowerBound = getMeanDist(traj2,selfTraj,numPts);
                    end
                    %Get pctError
                    pctError(index) = 100*(totalError-shiftError)/(totalError-lowerBound);
                    index = index + 1;
               end
           end
       end
    end
    
%% Plot results
    %Plot distributions
    figure; hold on
    histogram(pctError)

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