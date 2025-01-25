clear; clc; clf; close all;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = customRainbow;
    
%% Load Data 
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
    exclCh = [25 26 44 50 53 67 69 73 78 79 87 88 100 116 117 118 120 123];
    Data = getDataStruct20211210(Data,'getSpikes',true,'getSorts',false,'exclCh',exclCh,'exclZero',true,'getMarker',true,'centerMarker',false,'getKin',true);
    [Data] = cleanData20210706(Data);
    Data = Data([Data.trialStatus]==1);
    [Data,postureIDs] = labelPostures20210706(Data);
       
%% Bin data; smooth w 50ms std dev gaussian
    numTrials = size(Data,2);
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'GridReaching'};
    trialInclStates(1).inclStates = {{'state','Center Acquire','first',0},{'state','Last','last',0}};
    binWidth = 25;
    allFR = [];
    for trial = 1:numTrials
        trial
        [smoothFR,timestamps] = getStatesTraj20211210(Data(trial),trialInclStates,'smoothFR',binWidth,'timeRelToTrialStart',true);
        allFR = vertcat(allFR,smoothFR);
    end
        
 %% Collect all FR data; remove low FR channels; get z-score parameters; z-score data
    chMeans = mean(allFR);
    chStd = std(allFR);
    rmCh = find(chMeans < 1);
    
%% Get trajStruct 
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'GridReaching'};
    binWidth = 25;
    %Reach
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-100},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false,'chMeans',chMeans,'chStd',chStd,'rmCh',rmCh);   

%% Create distance distributions
    resultStruct = struct('cond1P',[],'cond1T',[],'cond2P',[],'cond2T',[],'dist',[]);
    numCond = size(trajStruct,2);
    numDraws = 100;
    numPts = 10;
    numSample = 4;
    structInd = 1;
    for cond1 = 1:numCond
        cond1
        for cond2 = cond1:numCond
            resultStruct(structInd).cond1P = trajStruct(cond1).posture;
            resultStruct(structInd).cond1T = trajStruct(cond1).target;
            resultStruct(structInd).cond2P = trajStruct(cond2).posture;
            resultStruct(structInd).cond2T = trajStruct(cond2).target;
            for i = 1:numDraws
                %Create traj1
                numTraj1 = size(trajStruct(cond1).allZSmoothFR,2);
                sampInd1 = randsample(numTraj1,numSample);
                traj1Struct = trajStruct(cond1).allZSmoothFR(sampInd1);
                traj1 = getAvgTraj20211210(traj1Struct,binWidth);
                %Create traj2
                if cond1 == cond2
                    numTrajRemaining = numTraj1 - numSample;
                    sampInd2 = randsample(numTrajRemaining,numSample);
                    remainingInd = setdiff(1:numTraj1,sampInd1);
                    sampInd2 = remainingInd(sampInd2);
                else
                    numTraj2 = size(trajStruct(cond2).allZSmoothFR,2);
                    sampInd2 = randsample(numTraj2,numSample);
                end
                traj2Struct = trajStruct(cond2).allZSmoothFR(sampInd2);
                traj2 = getAvgTraj20211210(traj2Struct,binWidth);
                %Get alpha
                alpha = getOptimalAlpha(traj1,traj2,numPts);
                %Shift traj2 by alpha
                traj2 = traj2 + alpha;
                %Get residual distance
                dist = getMeanDist(traj1,traj2,numPts);
                resultStruct(structInd).dist(i) = dist;
            end
            structInd = structInd + 1;
        end
    end
    
%% Collect lists of condition variables
cond1P = [resultStruct.cond1P]; cond1T = [resultStruct.cond1T]; cond2P = [resultStruct.cond2P]; cond2T = [resultStruct.cond2T];

%% Plot results
figure
condVec = [1 1 1 1];
dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
histogram(dist)

hold on

condVec = [2 2 2 2];
dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
histogram(dist)

legend('dist1','dist2')
xlabel('Mean Dist (Hz)')
ylabel('Count')

%% Check to see whether same condition distributions are overlapping 
for posture = 1:7
    for target = 1:8
        condVec = [posture target posture target];
        if any(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4))
            dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
            histogram(dist)
            hold on
        end
    end
end

xlabel('Mean Dist (Hz)')
ylabel('Count')

%% They are, so create chance dist out of them
chanceDist = [];
for posture = 1:7
    for target = 1:8
        condVec = [posture target posture target];
        if any(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4))
            dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
            chanceDist = [dist,chanceDist];
        end
    end
end

%% Compare conditions to chance dist

figure
histogram(chanceDist)
hold on

for target = 1:8
    for posture1 = 1:7
        for posture2 = posture1:7
            condVec = [posture1 target posture2 target];
            if any(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4))
                dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
                histogram(dist)
            end
        end
    end
end

%% Combine all test dist
testDist = [];
for target = 1:8
    for posture1 = 1:7
        for posture2 = posture1:7
            condVec = [posture1 target posture2 target];
            if any(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4))
                dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
                testDist = [dist,testDist];
            end
        end
    end
end

%% Plot chance and test
figure
histogram(chanceDist,'Normalization','probability')
hold on
histogram(testDist,'Normalization','probability')
% %% Test local functions
%     close all
%      
%     traj1 = [1 1;2 2];
% %     traj2 = traj1 + [1 1];
%     traj2 = [1 2; 2 4];
%     numPts = 2;
%     
%     figure
%     plot(traj1(:,1),traj1(:,2))
%     hold on
%     plot(traj2(:,1),traj2(:,2))
%     
%     getMeanDist(traj1,traj2,numPts)
%     alpha = getOptimalAlpha(traj1,traj2,numPts)
%     
%     figure
%     traj2 = traj2 + alpha;
%     plot(traj1(:,1),traj1(:,2))
%     hold on
%     plot(traj2(:,1),traj2(:,2))
%     
%     dist = getMeanDist(traj1,traj2,numPts)
%% Define local functions
    function alpha = getOptimalAlpha(traj1,traj2,numPts)
        numDims = size(traj1,2);
        alpha = NaN(1,numDims);
        for dim = 1:numDims
            alpha(dim) = (1/numPts)*sum(traj1(1:numPts,dim)-traj2(1:numPts,dim));
        end
    end
    
    function dist = getMeanDist(traj1,traj2,numPts)
        dist = 0;
        for i = 1:numPts
            dist = vecnorm(traj1(i,:)-traj2(i,:)) + dist;
        end
        dist = dist/numPts;
    end