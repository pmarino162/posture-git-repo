clear; clc; clf; close all;

%% Setup Save
    savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20210624\';
    saveFig = true;
    
%% Set parameters
    dateStr = '20210623';
    timeWindow = [100 500];
    getSorts = false;

%% Load Data 
    switch dateStr
        %6/22/21
        case '20210622'
            load('D:\Animals\Earl\2021\06\20210622\Earl20210622_Delayed Center Out_SI_SW_translated.mat');
            config = 'earlRHM1';
        %6/23/21
        case '20210623'
            load('D:\Animals\Earl\2021\06\20210623\Earl20210623_Delayed Center Out_SI_translated.mat');
            config = 'earlRHPMd';
    end
    
%% Get Data Struct 
    exclCh = [];
    Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'getMarker',true,'getKin',true);
    % Keep only successful trials 
    Data = Data([Data.trialStatus]==1);
    %If looking at delay tuning, keep only trials with appropriate delay
    %length
    kinData = [Data.kinData];
    delayLength = [kinData.delayLength];
    Data = Data(delayLength >= timeWindow(2));
    
%% Get Maps
    [tdt2adapterMap, adapter2cereportMap, cereport2arrayMap,array2anotomyMap,tdt2anatomyMap] = getMaps(config);
    
%% Create trajStruct
    condFields = {{'target','targetData','targetID'}};
    trajFields = {'allChannelSmoothedFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'Delayed Center Out 20210621'};
        trialInclStates(1).inclOccurrence = {'first'};
        trialInclStates(1).addTimeToBeginning = {timeWindow};
    %Delay
    trialInclStates(1).inclStates = {'Delay'}; 
    delayTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);    
    %Reach 
    trialInclStates(1).inclStates = {'Target Acquire'}; 
    reachTrajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates); 
    
%% Keep only trials from one posture 
%     trajStruct = trajStruct([trajStruct.reward]==2);

%% Keep only trials that fall within some time range

%% Get tuningData
    for period = {'delay','reach'}
        if strcmpi(period{1,1},'delay')
            trajStruct = delayTrajStruct;
        elseif strcmpi(period{1,1},'reach')
            trajStruct = reachTrajStruct;
        end
        % Get minimum number of trials satisfying criteria for each condition
        %Clear NANs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        minNumObs = zeros(1,size(trajStruct,2));
        for i = 1:size(trajStruct,2)
            if getSorts == true
                minNumObs(i) = size(trajStruct(i).allAllSortSmoothedFR,2);
            else
                minNumObs(i) = size(trajStruct(i).allAllChannelSmoothedFR,2);
            end
        end
        minNumObs = min(minNumObs);

        % Preallocate tuningData
        sortList = Data(1).spikes.sortList;
        numSorts = size(sortList,2);
        targetList = unique([trajStruct.target]);
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
            targetData = trajStruct([trajStruct.target]==target);
            for trial = 1:minNumObs
                if getSorts == true
                    trialFR = mean(targetData.allAllSortSmoothedFR(trial).traj);
                else
                    trialFR = mean(targetData.allAllChannelSmoothedFR(trial).traj);
                end
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
            for target = targetList
                startRow = minNumObs*(target-1)+1;
                endRow = startRow + minNumObs-1;
                y(startRow:endRow,1) = tuningData(sort).allData(:,target);
                x(startRow:endRow,:) = [ones(minNumObs,1),ones(minNumObs,1)*sind(targetAngles(target)),ones(minNumObs,1)*cosd(targetAngles(target))];
            end
            [B,bint,r,rint,stats] = regress(y,x);
            b0 = B(1,1); b1 = B(2,1); b2 = B(3,1); p = stats(3);
            tuningData(sort).MD = sqrt(b1.^2 + b2.^2);
            tuningData(sort).PD= atan2d(b1,b2);
            tuningData(sort).b0 = b0;
            tuningData(sort).p = p;
        end
        
        if strcmpi(period{1,1},'delay')
            delayTuningData = tuningData;
        elseif strcmpi(period{1,1},'reach')
            reachTuningData = tuningData;
        end
        clearvars tuningData trajStruct
    end
    
%% Plot results
    close all
    sigThreshold = 0.05;
    
    %PD Plot
    for period = {'delay','reach'}
        if strcmpi(period{1,1},'delay')
            tuningData = delayTuningData;
            color = [0 0 1];
        elseif strcmpi(period{1,1},'reach')
            tuningData = reachTuningData;
            color = [1 0 0];
        end
        figure('Position',[0 0 600 600])
        for sort = 1:numSorts
            PD = tuningData(sort).PD;
            MD = tuningData(sort).MD;
            quiver(0,0,MD*cosd(PD),MD*sind(PD),'LineWidth',2,'Color',color)
            hold on
        end
        sgtitle('Pref. Dirs scaled by Mod Depth')
        if saveFig 
            saveas(gcf,[savePath,period{1,1},'PDplot',dateStr,'.jpg'])
        end
    end

    %Tuning Curve Plot
    f = figure; f.Position = [10 10 700 800];
    [ha, pos] = tight_subplot(16,8,0,0,0);
    for period = {'delay','reach'}
        if strcmpi(period{1,1},'delay')
            tuningData = delayTuningData;
            color = [0 0 1];
        elseif strcmpi(period{1,1},'reach')
            tuningData = reachTuningData;
            color = [1 0 0];
        end
        for sort = 1:numSorts   
            channel = tuningData(sort).channel;
            row = tdt2anatomyMap{channel,2};
            col = tdt2anatomyMap{channel,3};
            axes(ha((row-1)*8 + col));
            avgFR = tuningData(sort).means;
            stdFR = tuningData(sort).SD;
            PD = tuningData(sort).PD;
            MD = tuningData(sort).MD;
            b0 = tuningData(sort).b0;
            p = tuningData(sort).p;
            Bfit = [b0;MD];
            x = [ones(numTargets,1),cosd(targetAngles-PD)];
            cosFit = x*Bfit;
            if p < sigThreshold
                plot(targetAngles,cosFit,'Color',color,'LineWidth',1.5)
            else
                plot(targetAngles,cosFit,'--','Color',color,'LineWidth',1.5)
            end
            hold on
            plot(targetAngles,avgFR,'.','Color',color,'MarkerSize',5);
            xticks([])
            yticks([])
            if strcmpi(period{1,1},'reach')
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
    end
%     numSigTuned = sum([tuningData.p] < sigThreshold);
%     sgtitle([num2str(numSigTuned),' Significantly Tuned'])
if saveFig 
    saveas(gcf,[savePath,'tuningCurves',dateStr,'.jpg'])
end

numChannels = 96;
delayMD = zeros(1,numChannels);
reachMD = zeros(1,numChannels);
channelList = [delayTuningData.channel];
for i = 1:numChannels
    if ismember(i,channelList)
        delayMD(i) = delayTuningData([delayTuningData.channel]==i).MD;
        reachMD(i) = reachTuningData([reachTuningData.channel]==i).MD;
    end
end

    createArrayImage(config,delayMD,dateStr)
    if saveFig 
        saveas(gcf,[savePath,'ArrayMapDelayMD',dateStr,'.jpg'])
    end

    createArrayImage(config,reachMD,dateStr)
    if saveFig 
        saveas(gcf,[savePath,'ArrayMapReachMD',dateStr,'.jpg'])
    end