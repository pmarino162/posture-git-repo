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
    
    trajStruct = reachTrajStruct;
    
%% Get Postures, Targets, and Channels
    postureList = unique([reachTrajStruct.posture]);
    targetList = unique([reachTrajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    numChannels = size(reachTrajStruct(1).avgAllChannelSmoothedFR.traj,2);
    
%% Do PCA on condition averages
    %Vertically concatenate trial averages
    allAvgs = [];
    for i = 1:size(trajStruct,2)
        traj = trajStruct(i).avgAllChannelSmoothedFR.traj;
        allAvgs = vertcat(allAvgs,traj);
    end
    [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);
    
%% Do LDA on all data 
    postureList = unique([trajStruct.posture]);
    targetList = unique([trajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    allTraj = []; allPostureLabels = []; allTargetLabels = [];
%     for i = 1:size(trajStruct,2)
%        target = trajStruct(i).target;
%        posture = trajStruct(i).posture;
%        numTraj = size(trajStruct(i).allAllChannelSmoothedFR,2);
%        for j = 1:numTraj
%             %Trajectory
%             traj = trajStruct(i).allAllChannelSmoothedFR(j).traj;
%             numSteps = size(traj,1);
%             allTraj = vertcat(allTraj,traj);
%             %Labels
%             trajTargetLabel = ones(numSteps,1).*target;
%             allTargetLabels = vertcat(allTargetLabels,trajTargetLabel);
%             trajPostureLabel = ones(numSteps,1).*posture;
%             allPostureLabels = vertcat(allPostureLabels,trajPostureLabel);
%         end
%     end

    
    for i = 1:size(trajStruct,2)
       target = trajStruct(i).target;
       posture = trajStruct(i).posture;
       numTraj = size(trajStruct(i).allAllChannelSmoothedFR,2);
       for j = 1:numTraj
            %Trajectory
            traj = trajStruct(i).allAllChannelSmoothedFR(j).traj;
            numSteps = size(traj,1);
            allTraj = vertcat(allTraj,traj(round(numSteps/2),:));
            %Labels
            trajTargetLabel = target;
            allTargetLabels = vertcat(allTargetLabels,trajTargetLabel);
            trajPostureLabel = posture;
            allPostureLabels = vertcat(allPostureLabels,trajPostureLabel);
        end
    end
    
    %Remove NaNs
    
    %Eliminate NaNs
    rmRows =[];
    for i = 1:size(allTraj,1)
       if any(isnan(allTraj(i,:)))
           rmRows = [rmRows,i];
       end
    end
    allTraj(rmRows,:) = []; allPostureLabels(rmRows,:) = []; allTargetLabels(rmRows,:);
    
    postureLDA = fisherLDA(allTraj, allPostureLabels);
    [postureLDA,~] = qr(postureLDA);
    postureLDA = postureLDA(:,1);
    
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    targetLDA = targetLDA(:,1:2);
%     postTargOrth = targetLDA;
    [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:3);
    
%% Add projections to trajStruct
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgPCA = (trajStruct(i).avgAllChannelSmoothedFR.traj-mu)*coeff;
       trajStruct(i).avgPostTargOrth = trajStruct(i).avgAllChannelSmoothedFR.traj*postTargOrth;
    end
    
%% Get Variance Explained in condition averages
    %Get total variance in condition averages
        totalVar = sum(var(allAvgs));
    %Get variance explained in each PCA dim
        pcaVarExpl = round((var((allAvgs-mu)*coeff)/totalVar)*100);
    %Get variance explained in each LDA dim
        ldaVarExpl = round((var(allAvgs*postTargOrth)/totalVar)*100);


%% Plot PSTHs
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

%% Preallocate posture tuning data 
    postureTuningData = struct('posture',[],'reachTuningData',[]);
    
%% Get tuningData
    postureList = unique([trajStruct.posture]);
    structInd = 1;
    for posture = postureList
        postureTuningData(structInd).posture = posture;
        tempTrajStruct = trajStruct([trajStruct.posture]==posture);
        % Get minimum number of trials satisfying criteria for each condition
        %Clear NANs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        minNumObs = zeros(1,size(tempTrajStruct,2));
        for i = 1:size(tempTrajStruct,2)
            minNumObs(i) = size(tempTrajStruct(i).allAllChannelSmoothedFR,2);
        end
        minNumObs = min(minNumObs);
        % Preallocate tuningData
        sortList = Data(1).spikes.sortList;
        numSorts = size(sortList,2);
        targetList = unique([tempTrajStruct.target]);
        numTargets = size(targetList,2);
        tuningData = struct('channel',nan,'sort',nan,'allData',nan(minNumObs,numTargets),'means',nan(1,numTargets),'SD',nan(1,numTargets),'MD',nan,'PD',nan,'b0',nan,'p',nan);
        tuningData = repmat(tuningData,numSorts,1);
        for i = 1:numSorts
            tuningData(i).channel = sortList(i).channel;
            tuningData(i).sort = sortList(i).sort;
        end
        % Get tuningData
        %allData
        targetInd = 1;
        for target = targetList
            targetData = tempTrajStruct([tempTrajStruct.target]==target);
            for trial = 1:minNumObs
                trialFR = mean(targetData.allAllChannelSmoothedFR(trial).traj(4:20,:));
                for sort = 1:numSorts
                    tuningData(sort).allData(trial,targetInd) = trialFR(sort);
                end
            end
            targetInd = targetInd + 1;
        end
        %mean and std
        for i = 1:size(tuningData,1)
           tuningData(i).means = nanmean(tuningData(i).allData);
           tuningData(i).SD = nanstd(tuningData(i).allData);
        end
        %Fit Tuning Curves
        targetAngles = transpose(45*(targetList-1)); 
        for sort = 1:numSorts
            y = nan(numTargets*minNumObs,1); x = nan(numTargets*minNumObs,3);
            for targetInd = 1:numTargets
                startRow = minNumObs*(targetInd-1)+1;
                endRow = startRow + minNumObs-1;
                y(startRow:endRow,1) = tuningData(sort).allData(:,targetInd);
                x(startRow:endRow,:) = [ones(minNumObs,1),ones(minNumObs,1)*sind(targetAngles(targetInd)),ones(minNumObs,1)*cosd(targetAngles(targetInd))];
            end
            [B,bint,r,rint,stats] = regress(y,x);
            b0 = B(1,1); b1 = B(2,1); b2 = B(3,1); p = stats(3);
            tuningData(sort).MD = sqrt(b1.^2 + b2.^2);
            PD = atan2d(b1,b2);
            if PD < 0
                PD = 360 - abs(PD);
            end
            tuningData(sort).PD= PD;
            tuningData(sort).b0 = b0;
            tuningData(sort).p = p;
        end
        postureTuningData(structInd).reachTuningData = tuningData;
        clearvars tuningData 
        structInd = structInd + 1;
    end
    
%% Plot PDs
    figure('Position',[0 0 600 600])
    color = [1 0 0];
    for neuron = 1:96
        PD = postureTuningData(1).reachTuningData(neuron).PD;
        MD = postureTuningData(1).reachTuningData(neuron).MD;
        Bv = [MD*cosd(PD),MD*sind(PD)];
        quiver(0,0,Bv(1),Bv(2),'LineWidth',2,'Color',color)
        hold on
    end
    sgtitle('Pref. Dirs scaled by Mod Depth')
%     saveas(gcf,[dirStr,'PDs\',profileMode,'.jpg'])
%     close 
    
%% Plot tuning curves
    close all
    sigThreshold = 0.05;
    postureListInd = 1;
    for posture1 = postureList
       for posture2 = postureList(postureListInd+1:end)
        f = figure; f.Position = [10 10 700 800];
        [ha, pos] = tight_subplot(12,8,0,0,0);
        postureInd = 1;
        for posture = [posture1,posture2]
            tuningData = postureTuningData([postureTuningData.posture]==posture).reachTuningData;
            %Get posture color
            if postureInd == 1
                color = [0 0 1];
            elseif postureInd == 2
                color = [1 0 0];
            end
            %Plot TC's
            for sort = 1:numSorts   
                row = floor(sort/8)+1;
                col = sort - 8*(row-1); 
                channel = tuningData(sort).channel;
                axes(ha((row-1)*8 + col));
                avgFR = tuningData(sort).means;
                stdFR = tuningData(sort).SD;
                PD = tuningData(sort).PD;
                MD = tuningData(sort).MD;
                b0 = tuningData(sort).b0;
                p = tuningData(sort).p;
                Bfit = [b0;MD];
                angleSpan = [0:45:315]';
                x = [ones(8,1),cosd(angleSpan-PD)];
                cosFit = x*Bfit;
                if p < sigThreshold
                    plot(angleSpan,cosFit,'Color',color,'LineWidth',1.5)
                else
                    plot(angleSpan,cosFit,'--','Color',color,'LineWidth',1.5)
                end
                hold on
                postureTargetList = [trajStruct([trajStruct.posture]==posture).target];
                targetAngles = (postureTargetList-1)*45;
                plot(targetAngles,avgFR,'.','Color',color,'MarkerSize',5);
                xticks([])
                yticks([])
                if postureInd == 2
                    yl = ylim;
                    xl = xlim;
                    scale = 0.7;
                    ylRange = yl(2)-yl(1);
                    ylMid = (yl(2)+yl(1))/2;
                    lowerLabel = round(ylMid-(ylRange/2)*scale);
                    upperLabel = round(ylMid+(ylRange/2)*scale);
                    text(0,lowerLabel,num2str(lowerLabel),'FontSize',8)
                    text(0,upperLabel,num2str(upperLabel),'FontSize',8)
                    text(xl(1)+(xl(2)-xl(1))*.7,upperLabel,['Ch' ,num2str(channel)],'FontSize',8)
                end
            end
            postureInd = postureInd + 1;
        end
       end
       postureListInd = postureListInd + 1;
    end
%     saveas(gcf,[dirStr,'TuningCurves\',profileMode,'.jpg'])
    close

%% Plot in PCA    
%     if size(postureList,2) == 2
%         cmap = cmap([1,5],:);
%     end
cmap = cmap([1,3,5],:);
    %Top 3 PCA
    figure
    xDim = 1; yDim =2; zDim = 3;
    postureInd = 1;
    for posture = postureList
        for target = targetList
            traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPCA;
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(traj(5,xDim),traj(5,yDim),traj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(traj(10,xDim),traj(10,yDim),traj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
        end
        postureInd = postureInd + 1;
    end
    xlabel(['PC ',num2str(xDim),' (',num2str(pcaVarExpl(xDim)),'%)']); ylabel(['PC ',num2str(yDim),' (',num2str(pcaVarExpl(yDim)),'%)']); zlabel(['PC ',num2str(zDim),' (',num2str(pcaVarExpl(zDim)),'%)']);
    grid on
    titleStr = ['Top 3 PC']; title(titleStr);
%     saveas(gcf,[dirStr,'PCA\',profileMode,'.jpg'])
%     saveas(gcf,[dirStr,'PCA\',profileMode,'.fig'])
%     close 
    
%% Plot in LDA 
    %Post Targ Orth - Actual
    figure
    xDim = 1; yDim =2; zDim = 3;
    postureInd = 1;
    for posture = postureList
        for target = targetList
            traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgPostTargOrth;
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(traj(5,xDim),traj(5,yDim),traj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(traj(10,xDim),traj(10,yDim),traj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
        end
        postureInd = postureInd + 1;
    end
    xlabel(['Posture LDA1 (',num2str(ldaVarExpl(xDim)),'%)']); ylabel(['Target LDA1 (',num2str(ldaVarExpl(yDim)),'%)']); zlabel(['Target LDA2 (',num2str(ldaVarExpl(zDim)),'%)']);
    grid on
    titleStr = ['LDA']; title(titleStr);
%     saveas(gcf,[dirStr,task,'ActualLDA.fig'])
%     saveas(gcf,[dirStr,'LDA\',profileMode,'.jpg'])
%     saveas(gcf,[dirStr,'LDA\',profileMode,'.fig'])
%     close 