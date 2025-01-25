clear; clc; clf; close all

%% Setup colormaps
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    
%% Generate Neural and Kinematic Data
    %Add poisson noise; make sure averages check out
    [kinStruct,neuronStruct,Data] = generateNeuralAndKinData;

%% Create whole trial, reaching, and holding traj structs
    %Parameters
    trajFields = {'allChannelSmoothedFR','marker'};
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
    for posture = [1,2]
       for target = 1:8
           if sum([trajStruct.posture]==posture & [trajStruct.target]==target)
               markerTraj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj;
               plot(markerTraj(:,1),markerTraj(:,2),'Color',tcmap(target,:));
               hold on
           end
       end
    end
    xlabel('x (mm)')
    ylabel('y (mm)')
    axis equal
    
%% Determine appropriate lag     
    resultStruct = struct('lag',[],'r2',[],'B',[],'p',[],'bOrth',[]);
    trajStruct = wholeTrajStruct;
    structInd = 1; 
    numCh = size(trajStruct(1).avgAllChannelSmoothedFR.traj,2);
    for lag = [-25:25:25]
        %Preallocate X and Y
        numObs = 0;
        numDims = 2;
        for i = 1:size(trajStruct,2)
            numTrials = size(trajStruct(i).allAllChannelSmoothedFR,2);
            for j = 1:numTrials
                neuralTraj = trajStruct(i).allAllChannelSmoothedFR(j).traj;
                numObs = numObs + size(neuralTraj,1);
            end
        end
        X = nan(numObs,3); Y = nan(numObs,numCh);
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
                X(k:k+numTrialObs-1,1) = ones(numTrialObs,1);
                X(k:k+numTrialObs-1,2:end) = markerInterp;
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
        
        %Perform cross-validated regression; save results
        numObs = size(X,1);
        ptn = cvpartition(numObs,'KFold',5);
        for fold = 1:5
            for channel = 1:numCh
                %Train
                trnIdx = training(ptn,fold);
                Xtrain = X(trnIdx,:);
                Ytrain = Y(trnIdx,channel);
                B = inv(Xtrain'*Xtrain)*Xtrain'*Ytrain;
                [b,bint,~,~,stats] = regress(Y(trnIdx,channel),X(trnIdx,:));
                %Test
                testIdx = test(ptn,fold);
                Xtest = X(testIdx,:);
                Ytest = Y(testIdx,channel);
                r2 = getR2(Ytest,Xtest*B);
                [F,p] = myFTest(Ytest,Xtest*B);
                %Save results
                resultStruct(structInd).r2(channel,fold) = r2;
                resultStruct(structInd).b(channel,:,fold) = b;
                resultStruct(structInd).p(channel,fold,:) = p;
                % qrb = qr(b(2:end,:));
                % resultStruct(structInd).bOrth(:,:,fold) = [qrb(:,1)./norm(qrb(:,1)),qrb(:,2)./norm(qrb(:,2))];
            end
        end
        resultStruct(structInd).lag = lag;
        structInd = structInd + 1;   
    end

% %% Plot r2 vs lag
% 
%     figure
%     subplot(2,1,1)
%         for i = 1:size(resultStruct,2)
%             lag = resultStruct(i).lag;
%             r2 = resultStruct(i).r2;
%             plot(lag,r2(:,1),'.','MarkerSize',5,'Color','b')
%             hold on
%         end
%         xlabel('lag (ms)')
%         ylabel('r2')
%     subplot(2,1,2)
%     for i = 1:size(resultStruct,2)
%         lag = resultStruct(i).lag;
%         r2 = resultStruct(i).r2;
%         plot(lag,r2(:,2),'.','MarkerSize',5,'Color','b')
%         hold on
%     end
%     xlabel('lag (ms)')
%     ylabel('r2')
%     

%% Choose best lag, rebuild X and Y for visualization
    structInd = 1; 
    numCh = size(trajStruct(1).avgAllChannelSmoothedFR.traj,2);
    lag = 0;
    %Preallocate X and Y
    numObs = 0;
    numDims = 2;
    for i = 1:size(trajStruct,2)
        numTrials = size(trajStruct(i).allAllChannelSmoothedFR,2);
        for j = 1:numTrials
            neuralTraj = trajStruct(i).allAllChannelSmoothedFR(j).traj;
            numObs = numObs + size(neuralTraj,1);
        end
    end
    X = nan(numObs,3); Y = nan(numObs,numCh);
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
            X(k:k+numTrialObs-1,1) = ones(numTrialObs,1);
            X(k:k+numTrialObs-1,2:end) = markerInterp;
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
    
%% Visualize fit for one neuron
    neuron = 1;
    B = resultStruct(2).b(neuron,1,:);
    B = reshape(B,1,3);
    lag = 0;
    figure
    plot3(X(:,2),X(:,3),Y(:,neuron),'.')
    ax = gca; xlim = ax.XLim; ylim = ax.YLim;
    [surfX,surfY] = meshgrid(xlim(1):10:xlim(2),ylim(1):10:ylim(2));
    surfZ = nan(size(surfX,1),size(surfX,2));
    for i = 1:size(surfX,1)
        for j = 1:size(surfX,2)
            surfZ(i,j) = B(1) + B(2)*surfX(i,j) + B(3)*surfY(i,j);
        end
    end
    hold on
    surf(surfX,surfY,surfZ);
    xlabel('Px (mm)'); ylabel('Py (mm)'); zlabel('FR (Hz)');
    
%% Decode arm position using PVA
B = resultStruct(2).b(:,:,1);
BMPinv = 
BMPinvMatlab = pinv(B)


%% Decode arm position using OLE 

