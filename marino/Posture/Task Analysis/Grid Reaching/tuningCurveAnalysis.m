clear; clc; clf; close all;

%% Setup Save
%     savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Task Analysis\Grid Reaching\';
    savePath = 'D:\Animals\Earl\2021\07\20210706\Daily Plots\Tuning Data\';
    saveFig = true;
    
%% Set parameters
    reachTimeWindow = [100 600];
    delayTimeWindow = [100 500];
    getSorts = false;
    
%% Load Data
    dateStr = '20210706';
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat')
%     load('D:\Animals\Earl\2021\07\20210702\Earl20210702_gridReaching_SI_translated.mat')
%     load('D:\Animals\Earl\2021\07\20210701\Earl20210701_gridReaching_SI_translated.mat')
%     load('D:\Animals\Earl\2021\06\20210630\Earl20210630_gridReaching_SI_translated.mat')
%     load('D:\Animals\Earl\2021\06\20210626\Earl20210626_gridReaching_SI_translated.mat');
%     load('D:\Animals\Earl\2021\06\20210628\Earl20210628_gridReaching_SI_translated.mat')
%     load('D:\Animals\Earl\2021\06\20210629\Earl20210629_gridReaching_SI_translated.mat')
%     config = 'earlRHM1';
    exclCh = [];

    if strcmpi(dateStr,'20210701')
        getSorts = false;
        exclCh = [1 6 24 33 35 36 40 43 44 69 52 53 57 55 58 67 69 70 71 77:80 82 87 88 91 95 96 119 120 121];        
    end

    
%% Get Data Struct 
    Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'exclZero',false,'getMarker',true,'getKin',true);
    % Keep only successful trials 
    Data = Data([Data.trialStatus]==1);
    
%% Clean Data for certain days
    if strcmpi(dateStr,'20210701')
        [Data] = cleanData20210701(Data);
    elseif strcmpi(dateStr,'20210706')
        [Data] = cleanData20210706(Data);
    end
    
%% Label Postures
    [Data,postureIDs] = labelPostures20210706(Data);
    
%% Get Maps
%     [tdt2adapterMap, adapter2cereportMap, cereport2arrayMap,array2anotomyMap,tdt2anatomyMap] = getMaps(config);
    
%% Keep trials appropriate for reach and for delay
    kinData = [Data.kinData];
    delayLength = [kinData.delayLength];
    delayData = Data(delayLength >= delayTimeWindow(2));
    reachData = Data;
    
%% Get Reach Trajectories
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    if getSorts
        trajFields = {'allSortSmoothedFR'};
    else
        trajFields = {'allChannelSmoothedFR'};
    end
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'GridReaching'};
        trialInclStates(1).inclOccurrence = {'first'};
        %trialInclStates(1).addTimeToEnd = {0};
    %Delay
    trialInclStates(1).inclStates = {'Delay'}; 
    trialInclStates(1).addTimeToBeginning = {delayTimeWindow};
    delayTrajStruct = getTrajStruct(delayData,condFields,trajFields,trialInclStates);    
    %Reach 
    trialInclStates(1).inclStates = {'Target Acquire'}; 
    trialInclStates(1).addTimeToBeginning = {reachTimeWindow};
    % Keep only trials with length w/in 2 stdDevs of mean length
    [reachData] = excludeLengths(reachData,trialInclStates);
    reachTrajStruct = getTrajStruct(reachData,condFields,trajFields,trialInclStates); 
    
%% Preallocate posture tuning data 
    postureTuningData = struct('posture',[],'delayTuningData',[],'reachTuningData',[]);

%% Get tuningData
    postureList = unique([reachTrajStruct.posture]);
    structInd = 1;
    for posture = postureList
        postureTuningData(structInd).posture = posture;
        for period = {'delay','reach'}
            if strcmpi(period{1,1},'delay')
                trajStruct = delayTrajStruct([delayTrajStruct.posture]==posture);
            elseif strcmpi(period{1,1},'reach')
                trajStruct = reachTrajStruct([reachTrajStruct.posture]==posture);
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

                    
            if strcmpi(period{1,1},'delay')
                postureTuningData(structInd).delayTuningData = tuningData;
            elseif strcmpi(period{1,1},'reach')
                postureTuningData(structInd).reachTuningData = tuningData;
            end
            clearvars tuningData trajStruct
        end
        structInd = structInd + 1;
    end
    
%% Get change in PD 
    delayChangePD = struct('posture1',[],'posture2',[],'changePD',[],'meanChangePD',[]);
    reachChangePD = struct('posture1',[],'posture2',[],'changePD',[],'meanChangePD',[]);

    sigThreshold = 0.05;
    for period = {'delay','reach'}
        structInd = 1;
        postureListInd = 1;
        for posture1 = postureList
           for posture2 = postureList(postureListInd:end)
                if strcmpi(period{1,1},'delay')
                    posture1TuningData = postureTuningData([postureTuningData.posture]==posture1).delayTuningData;
                    posture2TuningData = postureTuningData([postureTuningData.posture]==posture2).delayTuningData;
                elseif strcmpi(period{1,1},'reach')
                    posture1TuningData = postureTuningData([postureTuningData.posture]==posture1).reachTuningData;
                    posture2TuningData = postureTuningData([postureTuningData.posture]==posture2).reachTuningData;
                end
                sigTunedInBothPostures = find([posture1TuningData.p]<sigThreshold & [posture2TuningData.p]<sigThreshold);

                changePDInd = 1;
                for i = sigTunedInBothPostures
                   changePD(changePDInd).channel = posture1TuningData(i).channel;
                   changePD(changePDInd).sort = posture1TuningData(i).sort;
                   PD1 = posture1TuningData(i).PD;
                   PD2 = posture2TuningData(i).PD;
                   PDdiff1 = abs(PD1-PD2);
                   PDdiff2 = 360-PDdiff1;
                   changePD(changePDInd).changePD = min(PDdiff1,PDdiff2);
                   changePDInd = changePDInd + 1;
                end

                if strcmpi(period{1,1},'delay')
                    delayChangePD(structInd).posture1 = posture1;
                    delayChangePD(structInd).posture2 = posture2;
                    delayChangePD(structInd).changePD = changePD;
                    delayChangePD(structInd).meanChangePD = mean([changePD.changePD]);
                elseif strcmpi(period{1,1},'reach')
                    reachChangePD(structInd).posture1 = posture1;
                    reachChangePD(structInd).posture2 = posture2;
                    reachChangePD(structInd).changePD = changePD;
                    reachChangePD(structInd).meanChangePD = mean([changePD.changePD]);
                end
                structInd = structInd + 1;
                clearvars changePD
           end
           postureListInd = postureListInd + 1;
        end
    end

    edges = [0:5:360];
    for period = {'delay','reach'}
        if strcmpi(period{1,1},'delay')
            changePDStruct = delayChangePD;
        elseif strcmpi(period{1,1},'reach')
            changePDStruct = reachChangePD;
        end
        for i = 1:size(changePDStruct,2)
            figure
            posture1 = changePDStruct(i).posture1;
            posture2 = changePDStruct(i).posture2;
            histogram([changePDStruct(i).changePD.changePD],edges)
            meanChangePD = mean([changePDStruct(i).changePD.changePD]);
            hold on
            ylimits = ylim;
            line([meanChangePD, meanChangePD],ylimits,'Color','red','LineWidth',2)
            text(meanChangePD+3,ylimits(2)*.75,['mean = ',num2str(meanChangePD)])
            xlabel('Change in PD (degrees)')
            ylabel('Count')
            title([period,' changePDHist','Posture ',num2str(posture1),'vs ',num2str(posture2),' ',dateStr])
            if saveFig 
                saveas(gcf,[savePath,period{1,1},'ChangePDHist','Posture',num2str(posture1),'vs',num2str(posture2),'_',dateStr,'.jpg']);
            end
            close
        end
    end

%% Plot tuning curves
    close all
    sigThreshold = 0.05;
    for period = {'delay','reach'}
        postureListInd = 1;
        for posture1 = postureList
           for posture2 = postureList(postureListInd:end)
            f = figure; f.Position = [10 10 700 800];
            [ha, pos] = tight_subplot(16,8,0,0,0);
            postureInd = 1;
            for posture = [posture1,posture2]
                %Get tuning data
                if strcmpi(period{1,1},'delay')
                    tuningData = postureTuningData([postureTuningData.posture]==posture).delayTuningData;
                elseif strcmpi(period{1,1},'reach')
                    tuningData = postureTuningData([postureTuningData.posture]==posture).reachTuningData;
                end
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
        %             row = floor(channel/8)+1;
        %             col = channel - 8*(row-1);
        %             row = tdt2anatomyMap{channel,2};
        %             col = tdt2anatomyMap{channel,3};
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
        %             x = [ones(numTargets,1),cosd(targetAngles-PD)];
                    cosFit = x*Bfit;
                    if p < sigThreshold
                        plot(angleSpan,cosFit,'Color',color,'LineWidth',1.5)
                    else
                        plot(angleSpan,cosFit,'--','Color',color,'LineWidth',1.5)
                    end
                    hold on
                    postureTargetList = [reachTrajStruct([reachTrajStruct.posture]==posture).target];
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
            if saveFig 
                saveas(gcf,[savePath,period{1,1},'tuningCurves','Posture',num2str(posture1),'vs',num2str(posture2),'_',dateStr,'.jpg']);
            end
            close
           end
           postureListInd = postureListInd + 1;
        end
    end

%     numSigTuned = sum([tuningData.p] < sigThreshold);
%     sgtitle([num2str(numSigTuned),' Significantly Tuned'])

