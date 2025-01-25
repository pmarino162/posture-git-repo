clear; clc; clf; close all;

%% Set parameters 
    numCondTrials = 25;
    numNeurons = 96;
    ts = 25; %ms
    postureList = 1:2;
    targetList = 1:8;
    targetDist = 100; %mm
    reachTime = 500; %ms
    postureLocations = [0,0; 100,0];
    profileMode = 'bellVel';

%% Setup save dir   
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20211013 - moran model simulation\';

%% Create 1D position and velocity profiles
    time = 0:ts:reachTime;
    numSteps = size(time,2);
    
    switch profileMode
        case 'constVel'
            %Constant velocity
            velProf = (targetDist/reachTime)*ones(1,numSteps); %mm/ms
            posProf = zeros(1,numSteps);
            for i = 2:numSteps
               posProf(i) = posProf(i-1) + velProf(i-1)*ts; 
            end
        case 'bellVel'
            %Bell-shaped velocity profile
            a = pi/(targetDist^2);
            velProf = exp(-a*(time-150).^2); %mm/ms
            posProf = zeros(1,numSteps);
            for i = 2:numSteps
               posProf(i) = posProf(i-1) + velProf(i-1)*ts; 
            end
        case 'constVelStatPos'
            %Constant velocity; static position
            velProf = (targetDist/reachTime)*ones(1,numSteps); %mm/ms
            posProf = zeros(1,numSteps);
        case 'bellVelStatPos'
            %Bell-shaped velocity profile; static position
            a = pi/(targetDist^2);
            velProf = exp(-a*(time-150).^2); %mm/ms
            posProf = zeros(1,numSteps);
    end
    
    figure
    subplot(2,1,1)
        plot(time,velProf)
        xlabel('time (ms)')
        ylabel('speed (mm/ms)')
    subplot(2,1,2)
        plot(time,posProf)
        xlabel('time (ms)')
        ylabel('distance (mm)')
    close
        
%% Get position and velocity trajectories
    kinStruct = struct('posture',[],'target',[],'time',[],'position',[],'velocity',[]);
    structInd = 1;
    for posture = postureList
        for target = targetList
            theta = (target-1)*45;
            position = zeros(numSteps,2);
            velocity = zeros(numSteps,2);
            startPos = postureLocations(posture,:);   
            for i = 1:numSteps
                position(i,:) = startPos + [posProf(i)*cosd(theta), posProf(i)*sind(theta)];
                velocity(i,:) = [velProf(i)*cosd(theta), velProf(i)*sind(theta)];
            end
            kinStruct(structInd).posture = posture;
            kinStruct(structInd).target = target;
            kinStruct(structInd).time = time;
            kinStruct(structInd).position = position;
            kinStruct(structInd).velocity = velocity;
            structInd = structInd + 1;
        end
    end
    
    for posture = postureList
        for target = targetList
            time = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).time;
            position = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).position;
            velocity = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).velocity;
            
            figure
            subplot(2,2,1)
                plot(time,position(:,1)); xlabel('time (ms)'); ylabel('x (mm)')
            subplot(2,2,2)
                plot(time,position(:,2)); xlabel('time (ms)'); ylabel('y (mm)')
            subplot(2,2,3)
                plot(time,velocity(:,1)); xlabel('time (ms)'); ylabel('V_x (mm/ms)')
            subplot(2,2,4)
                plot(time,velocity(:,2)); xlabel('time (ms)'); ylabel('V_y (mm/ms)')
            sgtitle([profileMode,' T',num2str(target),' P',num2str(posture)])
%             saveas(gcf,[dirStr,'Kinematics\',profileMode,' T',num2str(target),' P',num2str(posture),' kinematics.jpg'])
            close 
        end
    end
    
%% Generate Neuron Properties
    neuronStruct = struct('Bv',[],'Bp',[],'b0',[]);
    meanMagBv = 87.5*(ts/1000);
    stdMagBv = meanMagBv/5;
    meanMagBp = .2*(ts/1000);
    stdMagBp = meanMagBp/5;
    meanb0 = 22.5*(ts/1000);
    stdb0 = 5*(ts/1000);
    for neuron = 1:numNeurons
       %Preferred vel direction
       PD = randn(1,2);
       PD = PD./norm(PD);
       %Preferred pos direction
       PP = randn(1,2);
       PP = PP./norm(PP);
       %Vel mod depth
       Mv = normrnd(meanMagBv,stdMagBv);
       %Pos mod depth
       Mp = normrnd(meanMagBp,stdMagBp);
       %Baseline FR
       b0 = normrnd(meanb0,stdb0); 
       %Model params
       neuronStruct(neuron).Bv = Mv.*PD;
       neuronStruct(neuron).Bp = Mp.*PP;
       neuronStruct(neuron).b0 = b0;
    end

% %% PDs
%     figure('Position',[0 0 600 600])
%     color = [1 0 0];
%     for neuron = 1:numNeurons
%         Bv = neuronStruct(neuron).Bv;
%         quiver(0,0,Bv(1),Bv(2),'LineWidth',2,'Color',color)
%         hold on
%     end
%     sgtitle('Pref. Dirs scaled by Mod Depth')

%% Generate Data Struct 
    stateData = struct('stateNames',{},'stateTransitions',[]);
    Data = struct('trialNum',[],'trialName','','trialStatus',[],'targetData',[],'conditionData',[],'stateData',stateData);
    sortList = struct('channel',[],'sort',[]);
    for i = 1:numNeurons
       sortList(i).channel = i;
       sortList(i).sort = 'all';
    end
    structInd = 1;
    trialNum = 0;
    for posture = postureList
        for target = targetList
            velocity = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).velocity;
            position = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).position;            
            for trial = 1:numCondTrials
                allChannelSpikeBins = zeros(size(time,2),numNeurons);
                for neuron = 1:numNeurons
                    Bv = neuronStruct(neuron).Bv;
                    Bp = neuronStruct(neuron).Bp;
                    b0 = neuronStruct(neuron).b0;
                    noise = normrnd(0,.1,numSteps,1);
                    allChannelSpikeBins(:,neuron) = [velocity*Bv' + position*Bp' + b0] + noise;
                end
                Data(structInd).trialNum = trialNum + 1;
                Data(structInd).trialName = 'Center Out';
                Data(structInd).conditionData.postureID = posture;
                Data(structInd).targetData.targetID = target;
                Data(structInd).stateData(1).stateNames = {'Reach','Success'};
                Data(structInd).stateData(1).stateTransitions = [1,2;0,reachTime];
                Data(structInd).spikes.binTimes = time;
                Data(structInd).spikes.allChannelSpikeBins = allChannelSpikeBins;
                Data(structInd).spikes.sortList = sortList;
                structInd = structInd + 1;
            end
        end
    end
    
%% Get trajStruct
    trajFields = {'allChannelSmoothedFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'posture','conditionData','postureID'},{'target','targetData','targetID'}};
    trialInclStates(1).trialName = {'Center Out'};
    trialInclStates(1).inclStates = {'Reach'};
    trialInclStates(1).inclOccurrence = {'first'};
    trialInclStates(1).addTimeToBeginning = {0};
    trialInclStates(1).addTimeToEnd = {0};   
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,'matchConditions',true);

%% Get Postures, Targets, and Channels
    postureList = unique([trajStruct.posture]);
    targetList = unique([trajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgAllChannelSmoothedFR.traj,2);
    
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

%% Setup colormaps
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;

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
                trialFR = mean(targetData.allAllChannelSmoothedFR(trial).traj(2:10,:));
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
    for neuron = 1:numNeurons
        PD = postureTuningData(1).reachTuningData(neuron).PD;
        MD = postureTuningData(1).reachTuningData(neuron).MD;
        Bv = [MD*cosd(PD),MD*sind(PD)];
        quiver(0,0,Bv(1),Bv(2),'LineWidth',2,'Color',color)
        hold on
    end
    sgtitle('Pref. Dirs scaled by Mod Depth')
%     saveas(gcf,[dirStr,'PDs\',profileMode,'.jpg'])
    close 
    
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
    if size(postureList,2) == 2
        cmap = cmap([1,5],:);
    end
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
    close 
    
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
    saveas(gcf,[dirStr,'LDA\',profileMode,'.jpg'])
    saveas(gcf,[dirStr,'LDA\',profileMode,'.fig'])
    close 