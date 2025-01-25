clear; clc; clf; close all

%% Setup colormaps
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    
%% Generate Neural and Kinematic Data
    posLag = 125;
    velLag = 75;
    [kinStruct,neuronStruct,Data] = generateNeuralAndKinData(posLag,velLag);

%% Create whole trial, reaching, and holding traj structs
    %Parameters
    trajFields = {'allChannelSmoothedFR','marker','markerVel'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'Center Out'};
    %Reaching
    trialInclStates(1).inclStates = {'Reach'};
    trialInclStates(1).inclOccurrence = {'first'};
    trialInclStates(1).addTimeToBeginning = {0};
    trialInclStates(1).addTimeToEnd = {0}; 
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
    %Whole Trial
    trialInclStates(1).inclStates = {'Reach','Hold'};
    trialInclStates(1).inclOccurrence = {'first','first'};
    trialInclStates(1).addTimeToBeginning = {0,0};
    trialInclStates(1).addTimeToEnd = {0,0}; 
    wholeTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
    
%% Get Postures, Targets, and Channels
    postureList = unique([trajStruct.posture]);
    targetList = unique([trajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgAllChannelSmoothedFR.traj,2);

%% Plot PSTHs
    trajStruct = wholeTrajStruct;
    f = figure; f.Position = [5 5 1400 700];
    [ha, pos] = tight_subplot(8,8,0.05,0.05,.05);
    channelInd = 1;
    for channel = 1:64
       axes(ha(channelInd));
       %Populate Figure
       for target = targetList
           for posture = 1
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgAllChannelSmoothedFR.traj;
                time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgAllChannelSmoothedFR.timestamps;
                hold on
                plot(time,traj(:,channel),'Color',tcmap(target,:),'LineWidth',2);
           end
       end
       channelInd = channelInd + 1;
       ax = gca;
       ax.TickDir = 'out';
       xticks([0 500])
       ylimits = ylim;
       if ylimits(1) > 0
           ax.YLim(1) = 0;
       end
       yticks([0 ylimits(2)])
       yticklabels({'0', [num2str(round(ylimits(2))),' Hz']})
       set(gca,'fontname','arial')
        xticklabels({'0','500 ms'})
    end
%     saveas(gcf,[dirStr,'PSTH\',profileMode,'.jpg'])

%% Determine appropriate lag     
    resultStruct = struct('posLag',[],'velLag',[],'r2px',[],'r2py',[],'r2vx',[],'r2vy',[],'B',[],'Bp',[],'Bv',[]);
    structInd = 1;
    numFolds = 5;
    numDims = 4;
    numCh = 96;
    for posLag = [0:25:200]
        posLag
        for velLag = [0:25:200]
            velLag
            %Preallocate
            numObs = 0;        
            for i = 1:size(trajStruct,2)
                numTrials = size(trajStruct(i).allAllChannelSmoothedFR,2);
                for j = 1:numTrials
                    neuralTraj = trajStruct(i).allAllChannelSmoothedFR(j).traj;
                    numObs = numObs + size(neuralTraj,1);
                end
            end
            X = nan(numObs,numDims); Y = nan(numObs,numCh);
            
            %Fill X and Y
            k = 1;
            for i = 1:size(trajStruct,2)
                %All Data
                numTrials = size(trajStruct(i).allAllChannelSmoothedFR,2);
                for j = 1:numTrials
                    %Get neural traj and time
                    neuralTraj = trajStruct(i).allAllChannelSmoothedFR(j).traj;
                    neuralTime = trajStruct(i).allAllChannelSmoothedFR(j).timestamps;
                    numTrialObs = size(neuralTime,2);
                    Y(k:k+numTrialObs-1,:) = neuralTraj;
                    %Get marker traj
                    markerTraj = trajStruct(i).allMarker(j).traj;
                    markerTime = trajStruct(i).allMarker(j).timestamps;
                    markerInterp = interp1(markerTime,markerTraj,neuralTime+posLag);
                    X(k:k+numTrialObs-1,1:2) = markerInterp;
                    markerVelTraj = trajStruct(i).allMarkerVel(j).traj;
                    markerVelTime = trajStruct(i).allMarkerVel(j).timestamps;
                    markerVelInterp = interp1(markerVelTime,markerVelTraj,neuralTime+velLag);
                    X(k:k+numTrialObs-1,3:4) = markerVelInterp;
                    %Update k
                    k = k + numTrialObs;
                end
            end
        
            %Eliminate NaNs
            rmRows =[];
            for i = 1:size(X,1)
               if any(isnan(X(i,:))) || any(isnan(Y(i,:)))
                   rmRows = [rmRows,i];
               end
            end
            X(rmRows,:) = []; Y(rmRows,:) = [];
            
            % Partition reaching data for CV
            numObs = size(X,1);
            ptn = cvpartition(numObs,'KFold',numFolds);
            %Add ones to X
            X = [ones(size(X,1),1),X];
            %Run model
            for fold = 1:numFolds
                trnIdx = training(ptn,fold);
                Xtrain = X(trnIdx,:);
                Ytrain = Y(trnIdx,:);
                %Train
                for channel = 1:numCh 
                    B = inv(Xtrain'*Xtrain)*Xtrain'*Ytrain(:,channel);
                    resultStruct(structInd).B(:,channel,fold) = B;
                end
                %Test
                testIdx = test(ptn,fold);
                testKin = X(testIdx,:);
                testNeural = Y(testIdx,:);
                B = resultStruct(structInd).B(:,:,fold);
                NKmap = pinv(B);
                predKin = testNeural*NKmap;
                resultStruct(structInd).r2px(1,fold) = getR2(testKin(:,2),predKin(:,2));
                resultStruct(structInd).r2py(1,fold) = getR2(testKin(:,3),predKin(:,3));
                resultStruct(structInd).r2vx(1,fold) = getR2(testKin(:,4),predKin(:,4));
                resultStruct(structInd).r2vy(1,fold) = getR2(testKin(:,5),predKin(:,5));
            end

            %Store results
            resultStruct(structInd).posLag = posLag;
            resultStruct(structInd).velLag = velLag;
            structInd = structInd + 1;
        end
    end
    