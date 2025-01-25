clear; clc; clf; close all;

%% Setup Save
%     savePath = 'D:\Animals\Earl\2021\07\20210728\Daily Plots\';
    saveFig = false;
   
%% Load Data
    dateStr = '20210805';
    getSorts = false;
    exclCh = [];
    load('D:\Animals\Earl\2021\08\20210805\17_N00_brainControl\Earl20210805_17_N00_brainControl_SI_translated.mat');
%      load('D:\Animals\Earl\2021\07\20210728\05_N00_brainControl\Earl20210728_05_N00_brainControl_SI_translated.mat')
%     load('D:\Animals\Earl\2021\07\20210713\03_brainControl\Earl20210713_03_brainControl_SI_translated.mat')
%% Get Data Struct 
    Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'exclZero',false,'getMarker',true,'getKin',true);
    % Keep only successful trials 
    Data = Data([Data.trialStatus]==1);
%     Data = Data(100:263);
%% Get Reach Trajectories
    condFields = {{'target','targetData','targetID'}};
    if getSorts
        trajFields = {'allSortSmoothedFR'};
    else
        trajFields = {'allChannelSmoothedFR'};
    end
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
%         
        
    
    trialInclStates(1).trialName = {'BCI Center Out'};   
    trialInclStates(1).inclStates = {'Step 1'}; 
    trialInclStates(1).inclOccurrence = {'first'};
    trialInclStates(1).addTimeToBeginning = {0};
    
%     trialInclStates(1).trialName = {'CenterOut_ForceBar_BC'}; 
%     trialInclStates(1).inclStates = {'BC Freeze','Touch Bar BC'}; 
%     trialInclStates(1).inclOccurrence = {'first','first'};
%     trialInclStates(1).addTimeToBeginning = {0,0};
    
    % Keep only trials with length w/in 2 stdDevs of mean length
    [Data] = excludeLengths(Data,trialInclStates);
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates); 
    
%% Get tuningData
    %Get min num Obs
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

 

%% PD plot 
    sigThreshold = 0.05;
    figure('Position',[0 0 600 600])
    for sort = 1:numSorts
        channel = tuningData(sort).channel;
        if ismember(channel,[33:96])
            color = [0 0 1];
        else 
            color = [1 0 0];
        end
        PD = tuningData(sort).PD;
        MD = tuningData(sort).MD;
        p = tuningData(sort).p;
        if p < sigThreshold
            quiver(0,0,MD*cosd(PD),MD*sind(PD),'LineWidth',2,'Color',color)
        else
        end
        hold on
    end
    xlim([-40 20])
    ylim([-25 20])
    sgtitle('Pref. Dirs scaled by Mod Depth')
    %Save
    if saveFig 
    saveas(gcf,[savePath,period{1,1},'PDplot',dateStr,'.jpg'])
    end


%% Plot tuning curves
    sigThreshold = 0.05;
    f = figure; f.Position = [10 10 700 800];
    [ha, pos] = tight_subplot(16,8,0,0,0);
    %Plot TC's
    for sort = 1:numSorts   
        row = floor(sort/8)+1;
        col = sort - 8*(row-1); 
        channel = tuningData(sort).channel;
        if ismember(channel,[33:96])
            color = [0 0 1];
        else 
            color = [1 0 0];
        end
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
        targetList = [trajStruct.target];
        targetAngles = (targetList-1)*45;
        plot(targetAngles,avgFR,'.','Color',color,'MarkerSize',5);
        xticks([])
        yticks([])

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

    %Save
    if saveFig 
        saveas(gcf,[savePath,period{1,1},'tuningCurves',dateStr,'.jpg']);
    end
%             close

%% Get number sig tuned 
    M1TuningData = tuningData([tuningData.channel]>=33 & [tuningData.channel]<=96);
    PMdTuningData = tuningData([tuningData.channel]<33 | [tuningData.channel]>96);
    
    numSigTunedM1 = sum([M1TuningData.p] < sigThreshold);
    numSigTunedPMd = sum([PMdTuningData.p] < sigThreshold);


