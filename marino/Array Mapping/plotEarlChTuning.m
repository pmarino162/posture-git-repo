clear; clc; clf; close all;

%% Setup Save
    savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20210624\';
    saveFig = true;
    
%% Set parameters
    startState = 'Target Acquire';
    startState = 'Delay';
    timeWindow = [100 500];
    getSorts = true;
    saveFig = false;
    savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Earl Array Analysis';

%% Get Data Struct 
    load('D:\Animals\Earl\2021\06\20210622\Earl20210622_Delayed Center Out_SI_SW_translated.mat');
    exclCh = [];
    Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'getMarker',true,'getKin',true);
    % Keep only successful trials 
    Data = Data([Data.trialStatus]==1);
    %If looking at delay tuning, keep only trials with appropriate delay
    %length
    kinData = [Data.kinData];
    delayLength = [kinData.delayLength];
    Data = Data(delayLength >= timeWindow(2));
    
%% Create trajStruct
    condFields = {{'target','targetData','targetID'}};
    trajFields = {'allSortSmoothedFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'Delayed Center Out 20210621'};
        trialInclStates(1).inclStates = {startState}; trialInclStates(1).inclOccurrence = {'first'};
        trialInclStates(1).addTimeToBeginning = {timeWindow};
%         trialInclStates(1).addTimeToEnd = {0};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);    
    
%% Keep only trials from one posture 
%     trajStruct = trajStruct([trajStruct.reward]==2);

%% Keep only trials that fall within some time range

%% Get minimum number of trials satisfying criteria for each condition
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

%% Preallocate tuningData
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
   
%% Get tuningData
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
        for targetInd = 1:numTargets
            startRow = minNumObs*(targetInd-1)+1;
            endRow = startRow + minNumObs-1;
            y(startRow:endRow,1) = tuningData(sort).allData(:,targetInd);
            x(startRow:endRow,:) = [ones(minNumObs,1),ones(minNumObs,1)*sind(targetAngles(targetInd)),ones(minNumObs,1)*cosd(targetAngles(targetInd))];
        end
        [B,bint,r,rint,stats] = regress(y,x);
        b0 = B(1,1); b1 = B(2,1); b2 = B(3,1); p = stats(3);
        tuningData(sort).MD = sqrt(b1.^2 + b2.^2);
        tuningData(sort).PD= atan2d(b1,b2);
        tuningData(sort).b0 = b0;
        tuningData(sort).p = p;
    end
        
%% Plot results
    sigThreshold = 0.05;
    %PD Plot
%     f(1)=figure('Position',[0 0 600 600])
%     f(1).Name = [fileName,'Vector']
    figure('Position',[0 0 600 600])
    for sort = 1:numSorts
        PD = tuningData(sort).PD;
        MD = tuningData(sort).MD;
        quiver(0,0,MD*cosd(PD),MD*sind(PD),'LineWidth',2,'Color','b')
        hold on
    end
    sgtitle('Pref. Dirs scaled by Mod Depth')
%     saveFigurePDF(f,'C:\Users\pmari\Documents\Energy Landscape\Presentations\7-5-19')
%     xlim([-3 3])
%     ylim([-3 3])
% 
    %Tuning Curve Plot
%     f(2)=figure
    numSortsInFig = 64;
    numSortsInRow = ceil(sqrt(numSortsInFig));
%     f(2).Name = [fileName,'Tuning']
    plotGroup = 0;
    for sort = 1:numSorts   
        sortGroup = floor((sort-1)/numSortsInFig)+1;
        if sortGroup > plotGroup
            plotGroup = sortGroup;
            figure('Position',[0 0 1200 800])
        end
        subplot(numSortsInRow,numSortsInRow,sort-(plotGroup-1)*numSortsInFig)
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
            plot(targetAngles,cosFit,'Color','r')
        else
            plot(targetAngles,cosFit,'--','Color','r')
        end
        hold on
        errorbar(targetAngles,avgFR,stdFR,'.k') % Necessary when there aren't 8 targets.
        plot(targetAngles,avgFR,'.','MarkerSize',2);
        xticks(targetAngles)
        xticklabels(num2str(targetAngles))
        xtickangle(45)
        channel = tuningData(sort).channel;
        sort = tuningData(sort).sort;
        title(['Ch ',num2str(channel),' Sort ',num2str(sort)])
    end

    numSigTuned = sum([tuningData.p] < sigThreshold);
    sgtitle([num2str(numSigTuned),' Significantly Tuned'])
    
%     if saveFig
%         if isempty(savePath)
%             savePath = uigetdir('Select Saving Location')
%         end
%         saveFigurePDF(f,savePath)
%     end
%     saveFigurePDF(f,'C:\Users\pmari\Documents\Energy Landscape\Presentations\7-5-19')


