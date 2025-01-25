clear; clc; clf; close all

%% Setup colormaps
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    
%% Generate Neural and Kinematic Data
    %Add poisson noise; make sure averages check out
    [kinStruct,neuronStruct,Data] = generateNeuralAndKinData(0,0);

%% Create whole trial, reaching, and holding traj structs
    %Parameters
    trajFields = {'allChannelSmoothedFR','marker','markerVel'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'Center Out'};
    %Whole Trial
    trialInclStates(1).inclStates = {'Reach','Hold'};
    trialInclStates(1).inclOccurrence = {'first','first'};
    trialInclStates(1).addTimeToBeginning = {0,0};
    trialInclStates(1).addTimeToEnd = {0,0}; 
    wholeTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
    %Reaching
    trialInclStates(1).inclStates = {'Reach'};
    trialInclStates(1).inclOccurrence = {'first'};
    trialInclStates(1).addTimeToBeginning = {0};
    trialInclStates(1).addTimeToEnd = {0}; 
    reachTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
    %Holding
    trialInclStates(1).inclStates = {'Hold'};
    trialInclStates(1).inclOccurrence = {'first'};
    trialInclStates(1).addTimeToBeginning = {0};
    trialInclStates(1).addTimeToEnd = {0}; 
    holdTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
    
%% Get Postures, Targets, and Channels
    postureList = unique([wholeTrajStruct.posture]);
    targetList = unique([wholeTrajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    numChannels = size(wholeTrajStruct(1).avgAllChannelSmoothedFR.traj,2);

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

%% Plot condition-averaged marker traces
    trajStruct = wholeTrajStruct;
    figure
    for posture = postureList
       for target = 1:8
           if sum([trajStruct.posture]==posture & [trajStruct.target]==target)
               markerTraj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj;
               plot(markerTraj(:,1),markerTraj(:,2),'Color',tcmap(target,:));
%                plot(markerTraj(:,1),markerTraj(:,2),'.','MarkerSize',10,'Color',tcmap(target,:));
               hold on
           end
       end
    end
    xlabel('x (mm)')
    ylabel('y (mm)')
    axis equal
    
%% Get X and Y for Hold and Reach
    lag = 0;
    numCh = size(trajStruct(1).avgAllChannelSmoothedFR.traj,2);
    for period = {'reach','hold'}
        %Get appropriate traj struct
        switch period{1,1}
            case 'reach'
                trajStruct = reachTrajStruct;
                numDims = 4;
            case 'hold'
                trajStruct = holdTrajStruct;
                numDims = 2;
        end
        %Preallocate X and Y
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
                markerInterp = interp1(markerTime,markerTraj,neuralTime+lag);
                X(k:k+numTrialObs-1,1:2) = markerInterp;
                if strcmp(period{1,1},'reach')
                    markerVelTraj = trajStruct(i).allMarkerVel(j).traj;
                    markerVelTime = trajStruct(i).allMarkerVel(j).timestamps;
                    markerVelInterp = interp1(markerVelTime,markerVelTraj,neuralTime+lag);
                    X(k:k+numTrialObs-1,3:4) = markerVelInterp;
                end
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
        
        %Save X and Y
        switch period{1,1}
            case 'reach'
                kinReach = X; neuralReach = Y;
            case 'hold'
                kinHold = X; neuralHold = Y;
        end
        clearvars X Y
    end

%% Partition reaching data for CV
    numFolds = 5;
    numObs = size(kinReach,1);
    ptn = cvpartition(numObs,'KFold',numFolds);
    
%% Create result struct
    resultStruct = struct('model','','r2px',[],'r2py',[],'r2vx',[],'r2vy',[],'B',[],'Bp',[],'Bv',[]);

%% Train and Test Simultaneous models
    structInd = 1;
    for model = {'PVA','PVA Pos Only','Reg','Reg Pos Only'}
        resultStruct(structInd).model = model{1,1};
        switch model{1,1}
            case 'PVA'
                X = kinReach; Y = neuralReach;
                X = [ones(size(X,1),1),X];
                %Perform cross-validated regression; save results
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
                
            case 'PVA Pos Only'
                X = kinReach(:,1:2); Y = neuralReach;
                X = [ones(size(X,1),1),X];
                %Perform cross-validated regression; save results
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
                end
                
            case 'Reg'
                X = neuralReach; Y = kinReach;
                X = [ones(size(X,1),1),X];
                %Perform cross-validated regression; save results
                for fold = 1:numFolds
                    %Train
                    trnIdx = training(ptn,fold);
                    Xtrain = X(trnIdx,:);
                    Ytrain = Y(trnIdx,:);
                    B = inv(Xtrain'*Xtrain)*Xtrain'*Ytrain;
                    resultStruct(structInd).B(:,:,fold) = B;
                    %Test
                    testIdx = test(ptn,fold);
                    Xtest = X(testIdx,:);
                    Ytest = Y(testIdx,:);
                    B = resultStruct(structInd).B(:,:,fold);
                    NKmap = B;
                    predKin = Xtest*NKmap;
                    resultStruct(structInd).r2px(1,fold) = getR2(Ytest(:,1),predKin(:,1));
                    resultStruct(structInd).r2py(1,fold) = getR2(Ytest(:,2),predKin(:,2));
                    resultStruct(structInd).r2vx(1,fold) = getR2(Ytest(:,3),predKin(:,3));
                    resultStruct(structInd).r2vy(1,fold) = getR2(Ytest(:,4),predKin(:,4));
                end
                
            case 'Reg Pos Only'
                X = neuralReach; Y = kinReach(:,1:2);
                X = [ones(size(X,1),1),X];
                %Perform cross-validated regression; save results
                for fold = 1:numFolds
                    %Train
                    trnIdx = training(ptn,fold);
                    Xtrain = X(trnIdx,:);
                    Ytrain = Y(trnIdx,:);
                    B = inv(Xtrain'*Xtrain)*Xtrain'*Ytrain;
                    resultStruct(structInd).B(:,:,fold) = B;
                    %Test
                    testIdx = test(ptn,fold);
                    Xtest = X(testIdx,:);
                    Ytest = Y(testIdx,:);
                    B = resultStruct(structInd).B(:,:,fold);
                    NKmap = B;
                    predKin = Xtest*NKmap;
                    resultStruct(structInd).r2px(1,fold) = getR2(Ytest(:,1),predKin(:,1));
                    resultStruct(structInd).r2py(1,fold) = getR2(Ytest(:,2),predKin(:,2));
                end
        end
        structInd = structInd + 1;   
    end

%% Train and test sequential models 
    for model = {'PVASeq','RegSeq'}
        resultStruct(structInd).model = model{1,1};
        switch model{1,1}
            case 'PVASeq'
                %Fit Bp on holding data
                X = kinHold; Y = neuralHold;
                X = [ones(size(X,1),1),X];
                for channel = 1:numCh 
                    Bp = inv(X'*X)*X'*Y(:,channel);
                    resultStruct(structInd).Bp(:,channel) = Bp;
                end
                %Remove position variability from neural signals
                X = kinReach; Y = neuralReach;
                X = [ones(size(X,1),1),X];
                kinReachPos = X(:,1:3);
                kinReachVel = X(:,[1,4,5]);
                Y = Y-kinReachPos*Bp;
                %Perform cross-validated regression; save results
                for fold = 1:numFolds
                    %Train
                    trnIdx = training(ptn,fold);
                    Xtrain = kinReachVel(trnIdx,:);
                    Ytrain = Y(trnIdx,:);
                    for channel = 1:numCh 
                        Bv = inv(Xtrain'*Xtrain)*Xtrain'*Ytrain(:,channel);
                        resultStruct(structInd).Bv(:,channel,fold) = Bv;
                    end
                    %Test
                    testIdx = test(ptn,fold);
                    testKin = kinReach(testIdx,:);
                    testNeural = neuralReach(testIdx,:);
                    Bp = resultStruct(structInd).Bp;
                    Bv = resultStruct(structInd).Bv(:,:,fold);
                    NPmap = pinv(Bp);
                    NVmap = pinv(Bv);
                    predPos = testNeural*NPmap;
                    testNeural = testNeural - [ones(size(testKin,1),1),testKin(:,1:2)]*Bp;
                    predVel = testNeural*NVmap;
                    resultStruct(structInd).r2px(1,fold) = getR2(testKin(:,1),predPos(:,2));
                    resultStruct(structInd).r2py(1,fold) = getR2(testKin(:,2),predPos(:,3));
                    resultStruct(structInd).r2vx(1,fold) = getR2(testKin(:,3),predVel(:,2));
                    resultStruct(structInd).r2vy(1,fold) = getR2(testKin(:,4),predVel(:,3));
                end
            case 'RegSeq'
                %Fit Bp on holding data
                X = neuralHold; Y = kinHold;
                X = [ones(size(X,1),1),X];
                Bp = inv(X'*X)*X'*Y;
                resultStruct(structInd).Bp = Bp;
                %Test
                X = neuralReach; Y = kinReach;
                X = [ones(size(X,1),1),X];
                for fold = 1:numFolds
                    testIdx = test(ptn,fold);
                    testNeural = X(testIdx,:);
                    testPos = Y(testIdx,1:2);
                    predKin = testNeural*Bp;
                    resultStruct(structInd).r2px(1,fold) = getR2(testPos(:,1),predKin(:,1));
                    resultStruct(structInd).r2py(1,fold) = getR2(testPos(:,2),predKin(:,2));
                end
        end
        structInd = structInd + 1;
    end
    
    
%% Visualize some model predictions
    model = resultStruct(1).model;
    B = resultStruct(1).B(:,:,fold);
    NKmap = pinv(B);

    %All targets
        posture = 1;
        %Position
        figure
        for target = 1:8
            neuralTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgAllChannelSmoothedFR.traj;
            markerTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgMarker.traj;
            predKin = neuralTraj*NKmap;
            plot(markerTraj(:,1),markerTraj(:,2),'b')
            hold on
            plot(predKin(:,2),predKin(:,3),'r')
        end
        xlabel('x (mm)')
        ylabel('y (mm)')
        legend('actual','pred')
        title(['Pos pred by ',model])
    
        %Velocity
        figure
        for target = 1:8
            neuralTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgAllChannelSmoothedFR.traj;
            markerVel = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgMarkerVel.traj;
            predKin = neuralTraj*NKmap;
            plot(markerVel(:,1),markerVel(:,2),'b')
            hold on
            plot(predKin(:,4),predKin(:,5),'r')
        end
        xlabel('Vx (mm/ms)')
        ylabel('Vy (mm/ms)')
        legend('actual','pred')
        title(['Vel pred by ',model])
    
    %One target vs time
    posture = 1;
    target = 8;
        %Position
        figure
        neuralTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgAllChannelSmoothedFR.traj;
        markerTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgMarker.traj;
        time = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgMarker.timestamps;
        predKin = neuralTraj*NKmap;
        sgtitle(['Pos pred by ',model])
        subplot(2,1,1)
            plot(time,markerTraj(:,1),'b')
            hold on
            plot(time,predKin(:,2),'r')
            xlabel('time (ms)')
            ylabel('x (mm)')
            legend('actual','pred')
        subplot(2,1,2)
            plot(time,markerTraj(:,2),'b')
            hold on
            plot(time,predKin(:,3),'r')
            xlabel('time (ms)')
            ylabel('y (mm)')
        
        %Velocity
        figure
        neuralTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgAllChannelSmoothedFR.traj;
        markerVel = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgMarkerVel.traj;
        time = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgMarkerVel.timestamps;
        predKin = neuralTraj*NKmap;
        sgtitle(['Vel pred by ',model])
        subplot(2,1,1)
            plot(time,markerVel(:,1),'b')
            hold on
            plot(time,predKin(:,4),'r')
            xlabel('time (ms)')
            ylabel('Vx (mm/ms)')
            legend('actual','pred')
        subplot(2,1,2)
            plot(time,markerVel(:,2),'b')
            hold on
            plot(time,predKin(:,5),'r')
            xlabel('time (ms)')
            ylabel('Vy (mm/ms)')
            
%% Plot Summary Prediction Performance 
    summaryPlotData = struct('model','','r2px',[],'r2py',[],'r2vx',[],'r2vy',[]);
    for i = 1:size(resultStruct,2)
       summaryPlotData(i).model = resultStruct(i).model;
       summaryPlotData(i).r2px = mean([resultStruct(i).r2px]);
       summaryPlotData(i).r2py = mean([resultStruct(i).r2py]);
    end
%Position Prediciton Performance 
    figure
    r2px = [summaryPlotData.r2px]';
    r2py = [summaryPlotData.r2py]';
    r2p = [r2px,r2py];
    modelName = {summaryPlotData.model};
    bar(r2p)
    title('Position Prediction')
    ylabel('r^2')

    xticklabels(modelName)
    xtickangle(45)
    legend('x','y')
    
    
%Velocity Prediction Performance 
    summaryPlotDataVel = struct('model','','r2vx',[],'r2vy',[]);
    j = 1;
    for i = [1,3,5]
       summaryPlotDataVel(j).model = resultStruct(i).model;
       summaryPlotDataVel(j).r2vx = mean([resultStruct(i).r2vx]);
       summaryPlotDataVel(j).r2vy = mean([resultStruct(i).r2vy]);
       j = j+1;
    end

    figure
    r2vx = [summaryPlotDataVel.r2vx]';
    r2vy = [summaryPlotDataVel.r2vy]';
    r2v = [r2vx,r2vy];
    modelName = {summaryPlotDataVel.model};
    bar(r2v)
    title('Velocity Prediction')
    ylabel('r^2')
    xticklabels(modelName)
    legend('x','y')