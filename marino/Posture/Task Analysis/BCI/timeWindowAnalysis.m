clear; clc; clf; close all

%% Load Data
%     load('D:\Animals\Earl\2021\07\20210729\05_N00_gradualTraining\Earl20210729_05_N00_gradualTraining_SI_translated.mat');
    load('D:\Animals\Earl\2021\08\20210802\04_N00_brainControl\Earl20210802_04_N00_brainControl_SI_translated.mat');
    Data = getDataStruct(Data,'getForce',true,'forceSetup','shoulder_posture_bci','getKin',true);

%% Get only successful data 
    sucData = Data([Data.trialStatus]==1);
    
%% Exlude trials that are too long
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'CenterOut_ForceBar_BC'};
        trialInclStates(1).trialName = {'CenterOut_ForceBar_BC'};
        trialInclStates(1).inclStates = {'BC Freeze','Touch Bar BC'}; trialInclStates(1).inclOccurrence = {'first','first'};
    [sucData] = excludeLengths(sucData,trialInclStates);
    
%% Get Reach Trajectories
    condFields = {{'target','targetData','targetID'}};
    trajFields = {'decoderGPFATraj','decoderCursorTraj','force'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'CenterOut_ForceBar_BC'};
        trialInclStates(1).inclStates = {'BC Freeze','Touch Bar BC'}; trialInclStates(1).inclOccurrence = {'first','first'};
    trajStruct = getTrajStruct(sucData,condFields,trajFields,trialInclStates);    

%% Separability by time bin
    numBins = 17;
    numObs = zeros(1,numBins);
    numFolds = 5;
    pctCorrectByBin = zeros(1,numBins);
    pctCorrectByBinCV = zeros(numBins,numFolds);
    
    for bin = 1:numBins
        %Preallocate
        obs = [];
        labels = [];
        %Get data
        for i = 1:size(trajStruct,2)
           targetID = trajStruct(i).target;
           allDecoderGPFATraj = trajStruct(i).allDecoderGPFATraj;
           for j = 1:size(allDecoderGPFATraj,2)
               traj = allDecoderGPFATraj(j).traj;
               numTimestamps = size(traj,1);
               if numTimestamps >= bin
                   obs = vertcat(obs,traj(bin,:));
                   labels = vertcat(labels,targetID);
               end
           end
        end
        %Classify
        [LDAModel,numClasses] = LDATrainModel(obs,labels);
        [predictedLabels,pctCorrect] = LDAClassify(obs,labels,numClasses,LDAModel);
        pctCorrectByBin(bin) = pctCorrect;
        numObs(bin) = size(obs,1);
        
        %Classify (Cross-Validated)
            numObs = size(obs,1);
            shufIdx = randperm(numObs);
            shufLabels = labels(shufIdx)';
            shufObs = obs(shufIdx,:);
            %Partition
            numObsInFold = floor(numObs/numFolds);
            shufLabels = shufLabels(1,1:numObsInFold*numFolds);
            shufObs = shufObs(1:numObsInFold*numFolds,:);
            for fold = 1:numFolds
                testInd = (fold-1)*numObsInFold+1:fold*numObsInFold;
                trainMask = ones(1,size(shufObs,1));
                trainMask(testInd) = 0;
                trainMask = logical(trainMask);
                trainObs = shufObs(trainMask,:);
                trainLabels = shufLabels(1,trainMask);
                testObs = shufObs(testInd,:);
                testLabels = shufLabels(1,testInd);
                %Fisher LDA
                [LDAModel,numClasses] = LDATrainModel(trainObs,trainLabels');
                [predictedLabels,pctCorrect] = LDAClassify(testObs,testLabels',numClasses,LDAModel);
                pctCorrectByBinCV(bin,fold) = pctCorrect;
            end
    end
    
    %Plot
    figure
    plot(pctCorrectByBin)
    xlabel('time bin (45ms bins)')
    ylabel('Classification Accuracy (%)')
    
    %Plot
    figure
    meanAccuracy = zeros(1,numBins);
    for bin = 1:numBins
       meanAccuracy(bin) = mean(pctCorrectByBinCV(bin,:));
       plot(ones(1,numFolds)*bin,pctCorrectByBinCV(bin,:),'.','MarkerSize',7,'Color','k')
       hold on
    end
    plot(1:numBins,meanAccuracy,'b','LineWidth',2)
    xlabel('time bin (45ms bins)')
    ylabel('Classification Accuracy (%)')
    
%% Distance from Baseline by time bin
    figure
    for target = 1:8
        traj = trajStruct([trajStruct.target]==target).avgDecoderGPFATraj.traj;
        dist = vecnorm(traj')';
        plot(dist)
        hold on
    end
    xlabel('time bin (45ms bins)')
    ylabel('Distance from Baseline (a.u.)')
    
%% Decoder cursor traj by bin 
    for target = 1:8
       traj = trajStruct([trajStruct.target]==target).avgDecoderCursorTraj.traj;
       figure
       subplot(2,1,1)
       plot(traj(:,1))
       ylabel('x (mm)')
       
       subplot(2,1,2)
       plot(traj(:,2))
       xlabel('time bin (45ms bins)')
       ylabel('y (mm)')
    end
    
%% Trajectories with time markers 
    figure
    hold on
    markerList = [5,7,9,11];
    markerColors = [0 0 0; 1 0 0; 1 0 1; 1 .5 0];
    for target = 1:8
       traj = trajStruct([trajStruct.target]==target).avgDecoderCursorTraj.traj;
       plot(traj(:,1),traj(:,2))
       markerInd = 1;
       for i = markerList
            plot(traj(i,1),traj(i,2),'.','MarkerSize',5,'Color',markerColors(markerInd,:));
            markerInd = markerInd + 1;
       end
    end
   xlabel('x (mm)')
   ylabel('y (mm)')