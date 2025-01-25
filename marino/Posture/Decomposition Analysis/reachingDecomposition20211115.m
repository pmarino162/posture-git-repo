clear; clc; clf; close all;
    
%% Setup save dir   
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20210930 - first model\';

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    
%% Load Data 
    binWidth = 25;
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
    exclCh =  [44 87 88 77 78 71 67 69 118];
    getSorts = false;
    Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'binWidth',binWidth,'exclZero',false,'getMarker',true,'centerMarker',false,'getKin',true);
    [Data] = cleanData20210706(Data);
    Data = Data([Data.trialStatus]==1);
    [Data,postureIDs] = labelPostures20210706(Data);
    allTrialPostures = [Data.conditionData];
    allTrialPostures = [allTrialPostures.postureID];
%          keepPostures = [1,2,4,5];
    keepPostures = [1,2,4,5,8,12,14];
    Data = Data(ismember(allTrialPostures,keepPostures));
    
%% Get Traj Struct   
    trajFields = {'allChannelSmoothedFR','marker'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'GridReaching'};
    trialInclStates(1).inclStates = {'Target Acquire'};
    trialInclStates(1).inclOccurrence = {'first'};
        trialInclStates(1).addTimeToBeginning = {0};
        trialInclStates(1).addTimeToEnd = {0}; 
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,'matchConditions',true);
        trialInclStates(1).addTimeToBeginning = {-250};
        trialInclStates(1).addTimeToEnd = {250}; 
    extTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,'matchConditions',true);


%% Plot condition-averaged marker traces
    for posture = [1,2,4,5,8,12,14]
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
    resultStruct = struct('lag',[],'r2',[],'b',[],'bOrth',[]);
    structInd = 1;
    for lag = [-200:25:200]
        lag
        
        %Preallocate
        numObs = 0;
        numCh = size(trajStruct(1).avgAllChannelSmoothedFR.traj,2);
        for i = 1:size(trajStruct,2)
            numTrials = size(trajStruct(i).allMarker,2);
            for j = 1:numTrials
                markerTraj = trajStruct(i).allMarker(j).traj(1:3:end,:);
                numObs = numObs + size(markerTraj,1);
            end
        end
        
        Y = nan(numObs,2); Ytime = nan(1,numObs); X = nan(numObs,numCh+1);
        k = 1;
        for i = 1:size(trajStruct,2)
            %All Data
            numTrials = size(trajStruct(i).allMarker,2);
            for j = 1:numTrials
                %Get marker traj; downsample; store
                markerTraj = trajStruct(i).allMarker(j).traj(1:3:end,:);
                markerTime = trajStruct(i).allMarker(j).timestamps(1:3:end);
                numTrialObs = size(markerTraj,1);
                Y(k:k+numTrialObs-1,:) = markerTraj;
                Ytime(1,k:k+numTrialObs-1) = markerTime;
%                 Y = vertcat(Y,markerTraj);
%                 Ytime = [Ytime,markerTime];
                %Get neural traj
                neuralTraj = extTrajStruct(i).allAllChannelSmoothedFR(j).traj;
                neuralTime = extTrajStruct(i).allAllChannelSmoothedFR(j).timestamps;
                neuralInterp = interp1(neuralTime,neuralTraj,markerTime+lag);
                X(k:k+numTrialObs-1,1) = ones(numTrialObs,1);
                X(k:k+numTrialObs-1,2:end) = neuralInterp;
%                 X = vertcat(X,neuralInterp);
                %Update k
                k = k + numTrialObs;
            end
        end
        numObs = size(X,1);
        %Perform cross-validated regression; save results
        ptn = cvpartition(numObs,'KFold',5);
        for fold = 1:5
            %Train
            trnIdx = training(ptn,fold);
            [bX,bintX,~,~,statsX] = regress(Y(trnIdx,1),X(trnIdx,:));
            [bY,bintY,~,~,statsY] = regress(Y(trnIdx,2),X(trnIdx,:));
            b = [bX,bY];
            %Test
            testIdx = test(ptn,fold);
            rX = getR2(Y(testIdx,1),X(testIdx,:)*bX);
            rY = getR2(Y(testIdx,2),X(testIdx,:)*bY);
            %Save results
            resultStruct(structInd).r2(fold,:) = [rX,rY];
            resultStruct(structInd).b(:,:,fold) = b;
            qrb = qr(b(2:end,:));
            resultStruct(structInd).bOrth(:,:,fold) = [qrb(:,1)./norm(qrb(:,1)),qrb(:,2)./norm(qrb(:,2))];
        end
        resultStruct(structInd).lag = lag;
        structInd = structInd + 1;   
    end

%% Plot r2 vs lag

    figure
    subplot(2,1,1)
        for i = 1:size(resultStruct,2)
            lag = resultStruct(i).lag;
            r2 = resultStruct(i).r2;
            plot(lag,r2(:,1),'.','MarkerSize',5,'Color','b')
            hold on
        end
        xlabel('lag (ms)')
        ylabel('r2')
    subplot(2,1,2)
    for i = 1:size(resultStruct,2)
        lag = resultStruct(i).lag;
        r2 = resultStruct(i).r2;
        plot(lag,r2(:,2),'.','MarkerSize',5,'Color','b')
        hold on
    end
    xlabel('lag (ms)')
    ylabel('r2')
    
%% Choose lag, Visualize some data and fits
    %Choose lag
    lag = -50;
    b = resultStruct(find([resultStruct.lag]==lag)).b;
    b = mean(b,3);
    %Build X and Y

    %Preallocate
    numObs = 0;
    numCh = size(trajStruct(1).avgAllChannelSmoothedFR.traj,2);
    for i = 1:size(trajStruct,2)
        numTrials = size(trajStruct(i).allMarker,2);
        for j = 1:numTrials
            markerTraj = trajStruct(i).allMarker(j).traj(1:3:end,:);
            numObs = numObs + size(markerTraj,1);
        end
    end
        
    Y = nan(numObs,2); Ytime = nan(1,numObs); X = nan(numObs,numCh+1);
    k = 1;
    for i = 1:size(trajStruct,2)
        %All Data
        numTrials = size(trajStruct(i).allMarker,2);
        for j = 1:numTrials
            %Get marker traj; downsample; store
            markerTraj = trajStruct(i).allMarker(j).traj(1:3:end,:);
            markerTime = trajStruct(i).allMarker(j).timestamps(1:3:end);
            numTrialObs = size(markerTraj,1);
            Y(k:k+numTrialObs-1,:) = markerTraj;
            Ytime(1,k:k+numTrialObs-1) = markerTime;
            %Get neural traj
            neuralTraj = extTrajStruct(i).allAllChannelSmoothedFR(j).traj;
            neuralTime = extTrajStruct(i).allAllChannelSmoothedFR(j).timestamps;
            neuralInterp = interp1(neuralTime,neuralTraj,markerTime+lag);
            X(k:k+numTrialObs-1,1) = ones(numTrialObs,1);
            X(k:k+numTrialObs-1,2:end) = neuralInterp;
            %Update k
            k = k + numTrialObs;
        end
    end
    
    
% %     Y = []; Ytime = []; X = [];
% %     for i = 1:size(trajStruct,2)
% %         %Get marker traj; downsample; store
% %         markerTraj = trajStruct(i).avgMarker.traj(1:3:end,:);
% %         markerTime = trajStruct(i).avgMarker.timestamps(1:3:end);
% %         Y = vertcat(Y,markerTraj);
% %         Ytime = [Ytime,markerTime];
% %         %Get neural traj
% %         neuralTraj = extTrajStruct(i).avgAllChannelSmoothedFR.traj;
% %         neuralTime = extTrajStruct(i).avgAllChannelSmoothedFR.timestamps;
% %         neuralInterp = interp1(neuralTime,neuralTraj,markerTime+lag);
% %         X = vertcat(X,neuralInterp);
% %     end
% %     numObs = size(X,1);
% %     X = [ones(numObs,1),X];

    %Single neuron visualization
    neuron = 10;
    
    %Model fit visualization
    figure
        plot(Y(:,1),X*b(:,1),'.')
        title('X position')
        xlabel('x (mm)')
        ylabel('x_{Pred} (mm)')
    figure
        plot(Y(:,2),X*b(:,2),'.')
        title('Y position')
        xlabel('y (mm)')
        ylabel('y_{Pred} (mm)')

%% Plot kinematic and predicted kinematic traces 
    lag = -50;
    for target = 1:8
        posture = 1;
    
        markerTraj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj;
        markerTime = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps;
%         markerTime = markerTime + lag;
    
        neuralTraj = extTrajStruct([extTrajStruct.posture]==posture & [extTrajStruct.target]==target).avgAllChannelSmoothedFR.traj;       
        neuralTime = extTrajStruct([extTrajStruct.posture]==posture & [extTrajStruct.target]==target).avgAllChannelSmoothedFR.timestamps;
        neuralInterp = interp1(neuralTime,neuralTraj,markerTime+lag);
        neuralTraj = [ones(size(neuralInterp,1),1),neuralInterp];
        
        
        markerPred = neuralTraj*b;
    
        figure
        sgtitle(['Target ',num2str(target)])
        subplot(2,1,1)
            plot(markerTime,markerTraj(:,1))
            hold on
            plot(markerTime,markerPred(:,1))
            xlabel('time (ms)')
            ylabel('x (mm)')
        subplot(2,1,2)
            plot(markerTime,markerTraj(:,2))
            hold on
            plot(markerTime,markerPred(:,2))
            xlabel('time (ms)')
            ylabel('y (mm)')
        saveas(gcf,['C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20211118 - posture regression during reaching 1\','Target ',num2str(target),' Prediction.jpg'])
        close all
    end
        
%% Choose coefficients from model with best cross-validated r2

%% Remove chosen dimensions from data and run PCA on null Space
%     qrb = qr(b(2:end,:));
%     bOrth = [qrb(:,1)./norm(qrb(:,1)),qrb(:,2)./norm(qrb(:,2))];
    
    
    bOrth = [b(2:end,1)/norm(b(2:end,1)),b(2:end,2)/norm(b(2:end,2))];
    postureNullBasis = null(bOrth');
    
    
    %Vertically concatenate trial averages
    allAvgs = [];
    for i = 1:size(trajStruct,2)
        traj = trajStruct(i).avgAllChannelSmoothedFR.traj;
        allAvgs = vertcat(allAvgs,traj);
    end
    nullProj = allAvgs*postureNullBasis;
    %Create null Basis
    [postureNullBasisPC,score,latent,tsquared,explained,mu] = pca(nullProj);
    postureNullBasis = postureNullBasis*postureNullBasisPC;
           
    
%% Add Projections to trajStruct and modelStruct
    for i = 1:size(trajStruct,2)
       trajStruct(i).postureAxes = (trajStruct(i).avgAllChannelSmoothedFR.traj)*bOrth;
       trajStruct(i).avgPCA = (trajStruct(i).avgAllChannelSmoothedFR.traj)*postureNullBasis;
    end    
    
%% Get Variance Explained in condition averages
    %Get total variance in condition averages
        totalVar = sum(var(allAvgs));

    %Get variance explained in each PCA dim
        postureVarExpl = (var(allAvgs*bOrth)/totalVar)*100;
        pcaVarExpl = (var(allAvgs*postureNullBasis)/totalVar)*100;
        varExpl = [postureVarExpl,pcaVarExpl];
        

    
    
%% PC Timecourses - multiple postures for each target

%     postureList = [1:5];
    postureList = [2,12];
%     postureList = [4,14];
%     postureList = [1,2,4,5];
    
   cmap = customSummer;
    if size(postureList,2) == 2
        cmap = cmap([1,5],:);
    end
%     
    for target = 1:8
            f = figure; f.Position = [10 10 1500 700];
                    sgtitle(['Target ',num2str(target)])
    for posture = postureList
        postureTraj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).postureAxes;
        nullTraj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPCA;
        traj = horzcat(postureTraj,nullTraj);
        time = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgAllChannelSmoothedFR.timestamps;
        for plotInd = 1:10
            subplot(2,5,plotInd)
            plot(time,traj(:,plotInd),'Color',cmap([find(postureList==posture)],:),'LineWidth',2);
            hold on
            if posture == postureList(end)
               xlabel('time (ms)')
               if plotInd < 3
                   ylabel(['Posture Regression ','(',num2str(varExpl(plotInd)),'%)'])
               else
                   ylabel(['Null Dim ',num2str(plotInd-2),' (',num2str(varExpl(plotInd)),'%)'])
               end  
            end
        end
    end
    
    %Get axis ranges
    for plotInd = 1:10
       ax = subplot(2,5,plotInd);
       range(plotInd) = ax.YLim(1,2)-ax.YLim(1,1);
       center(plotInd) = (ax.YLim(1,2)+ax.YLim(1,1))/2;
    end
    maxRange = max(range);
    for plotInd = 1:10
       ax = subplot(2,5,plotInd);
       ax.YLim = [center(plotInd)-maxRange/2 center(plotInd)+maxRange/2];
    end
    
    saveas(gcf,['C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20211118 - posture regression during reaching 1\','Target ',num2str(target),' NeuralTraj.jpg'])
    close all
    
    end
    
%% PC Timecourses - multiple targets and postures
    f = figure; f.Position = [10 10 1500 700];
%     postureList = [1,5];
    targetList = [3,7];
    postureList = [1,2,4,5];
    for posture = postureList
        for target = targetList
            postureTraj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).postureAxis;
            nullTraj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPCA;
            traj = horzcat(postureTraj,nullTraj);
            time = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgAllChannelSmoothedFR.timestamps;
            for plotInd = 1:10
                subplot(2,5,plotInd)
                if target == targetList(end)
%                     plot(time,traj(:,plotInd),'--','Color',cmap([find(postureList==posture)],:),'LineWidth',2);
                    plot(time,traj(:,plotInd),'--','Color',cmap(posture,:),'LineWidth',2);
                    hold on
                else
%                     plot(time,traj(:,plotInd),'Color',cmap([find(postureList==posture)],:),'LineWidth',2);
                    plot(time,traj(:,plotInd),'Color',cmap(posture,:),'LineWidth',2);
                    hold on
                end
                if posture == postureList(end)
                   xlabel('time (ms)')
                   if plotInd == 1
                       ylabel(['Posture LDA ','(',num2str(varExpl(plotInd)),'%)'])
                   else
                       ylabel(['Null Dim ',num2str(plotInd-1),' (',num2str(varExpl(plotInd)),'%)'])
                   end  
                end
            end
        end
    end
    
    %Get axis ranges
    for plotInd = 1:10
       ax = subplot(2,5,plotInd);
       range(plotInd) = ax.YLim(1,2)-ax.YLim(1,1);
       center(plotInd) = (ax.YLim(1,2)+ax.YLim(1,1))/2;
    end
    maxRange = max(range);
    for plotInd = 1:10
       ax = subplot(2,5,plotInd);
       ax.YLim = [center(plotInd)-maxRange/2 center(plotInd)+maxRange/2];
    end 

    
%             %Condition Averages
%                 %Get marker traj; downsample; store
%                 markerTraj = trajStruct(i).avgMarker.traj(1:3:end,:);
%                 markerTime = trajStruct(i).avgMarker.timestamps(1:3:end);
%                 Y = vertcat(Y,markerTraj);
%                 Ytime = [Ytime,markerTime];
%                 %Get neural traj
%                 neuralTraj = extTrajStruct(i).avgAllChannelSmoothedFR.traj;
%                 neuralTime = extTrajStruct(i).avgAllChannelSmoothedFR.timestamps;
%                 neuralInterp = interp1(neuralTime,neuralTraj,markerTime+lag);
%                 X = vertcat(X,neuralInterp);