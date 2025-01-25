clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = true;
    saveDir = 'D:\Posture Paper\20221123\Figures\Figure 1';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%% Load dataset
    %dataset = 'R20200221';
    dataset = 'E20210706';
    [Data,zScoreParams] = loadData(dataset);
    
%% Get trajStruct
    binWidth = 50;
    kernelStdDev = 50;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    switch dataset
        case 'E20210706'
             trialInclStates(1).trialName = {'GridReaching'};
             condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
             trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',200}};
        case 'R20200221'
             trialInclStates(1).trialName = {'Rocky Dissociation'};
             condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
             trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',200}};  
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    
    if strcmpi(dataset,'E20210706')
        trajStruct = trajStruct(ismember([trajStruct.posture],1:3));
    end
    
    
%% Get posture tuning data 
    postureTuningData = struct('posture',[],'tuningData',[]);
    postureList = unique([trajStruct.posture]);
    targetList = unique([trajStruct.target]);
    numTargets = size(targetList,2);
    numCh = size(Data(1).spikes,2);
    numTrials = size(Data,2);
    
    structInd = 1;
    for posture = postureList
        % Preallocate tuningData
        tuningData = struct('channel',nan,'allData',nan(numTrials,numTargets),'means',nan(1,numTargets),'SD',nan(1,numTargets),'MD',nan,'PD',nan,'b0',nan,'p',nan);
        tuningData = repmat(tuningData,numCh,1);
        for i = 1:numCh
            tuningData(i).channel = i;
        end
        % Get tuningData
        tempTrajStruct = trajStruct([trajStruct.posture]==posture);
        targetInd = 1;
        for target = targetList
            targetData = tempTrajStruct([tempTrajStruct.target]==target);
            for trial = 1:numel(targetData.allSmoothFR)
                trialFR = mean(targetData.allSmoothFR(trial).traj);
                for channel = 1:numCh
                    tuningData(channel).allData(trial,targetInd) = trialFR(channel);
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
        for channel = 1:numCh
            y = nan(numTargets,1); x = nan(numTargets,3);
            for targetInd = 1:numTargets
                y(targetInd,1) = tuningData(channel).means(1,targetInd);
                x(targetInd,:) = [1,sind(targetAngles(targetInd)),cosd(targetAngles(targetInd))];
            end
            [B,bint,r,rint,stats] = regress(y,x);
            b0 = B(1,1); b1 = B(2,1); b2 = B(3,1); p = stats(3);
            tuningData(channel).MD = sqrt(b1.^2 + b2.^2);
            PD = atan2d(b1,b2);
            if PD < 0
                PD = 360 - abs(PD);
            end
            tuningData(channel).PD= PD;
            tuningData(channel).b0 = b0;
            tuningData(channel).p = p;
        end
        %Fill Struct
        postureTuningData(structInd).posture = posture;
        postureTuningData(structInd).tuningData = tuningData;
        structInd = structInd + 1;
    end
    
%% Plot tuning curves
    sigThreshold = 0.10;
    numChannelGroups = ceil(numCh/64);
    for channelGroup = 1:numChannelGroups
        channelList = (channelGroup-1)*64+1:(channelGroup)*64;
        channelList = channelList(channelList <= numCh);
        f = figure; f.Position = [10 10 700 650];
        [ha, pos] = tight_subplot(8,8,0,0,0);
        postureInd = 1;
        postureListInd = 1;
        targetAngles = (targetList-1)*45;
        for posture = postureList
            %Get tuning data
            tuningData = postureTuningData([postureTuningData.posture]==posture).tuningData;
            %Plot TC's
            for channel = channelList
                row = floor(channel-(channelList(1)-1)/8)+1;
                col = channel-(channelList(1)-1) - 8*(row-1); 
                axes(ha((row-1)*8 + col));
                avgFR = tuningData(channel).means;
                stdFR = tuningData(channel).SD;
                PD = tuningData(channel).PD;
                MD = tuningData(channel).MD;
                b0 = tuningData(channel).b0;
                p = tuningData(channel).p;
                Bfit = [b0;MD];
                angleSpan = [0:45:315]';
                x = [ones(8,1),cosd(angleSpan-PD)];
                cosFit = x*Bfit;
                if p < sigThreshold
                    plot(angleSpan,cosFit,'Color',pcmap(posture,:),'LineWidth',1.5)
                else
                    plot(angleSpan,cosFit,'--','Color',pcmap(posture,:),'LineWidth',1.5)
                end
                hold on;
                plot(targetAngles,avgFR,'.','Color',pcmap(posture,:),'MarkerSize',5);
                xticks([]); yticks([])        
                if posture == postureList(end)
                    yl = ylim;
                    xl = xlim;
                    scale = 0.7;
                    ylRange = yl(2)-yl(1);
                    ylMid = (yl(2)+yl(1))/2;
                    %lowerLabel = round(ylMid-(ylRange/2)*scale);
                    upperLabel = ylMid+(ylRange/2)*scale;
                    text(xl(1),upperLabel,num2str(channel),'FontSize',8)
    %                 text(0,lowerLabel,num2str(lowerLabel),'FontSize',8)
    %                 text(0,upperLabel,num2str(upperLabel),'FontSize',8)
    %                 text(xl(1)+(xl(2)-xl(1))*.7,upperLabel,['Ch' ,num2str(channel)],'FontSize',8)
                end         
            end
            postureInd = postureInd + 1;
            postureListInd = postureListInd + 1;     
        end
        
        
    end
    
  

%% Plot tuning curves for select channels
    plotChannelList = [103,32,108,80];
    numPlotCh = numel(plotChannelList);
    sigThreshold = 0.10;
    f = figure; f.Position = [10 10 500 100];
    [ha, pos] = tight_subplot(1,numPlotCh,0,0,0);
    postureInd = 1;
    postureListInd = 1;
    
    targetAngles = (targetList-1)*45;

    for posture = postureList
        %Get tuning data
        tuningData = postureTuningData([postureTuningData.posture]==posture).tuningData;
        %Plot TC's
        channelInd = 1;
        for channel = plotChannelList
            row = 1;
            col = channelInd; 
            axes(ha(col));
            avgFR = tuningData(channel).means;
            stdFR = tuningData(channel).SD;
            PD = tuningData(channel).PD;
            MD = tuningData(channel).MD;
            b0 = tuningData(channel).b0;
            p = tuningData(channel).p;
            Bfit = [b0;MD];
            angleSpan = [0:45:315]';
            x = [ones(8,1),cosd(angleSpan-PD)];
            cosFit = x*Bfit;
            if p < sigThreshold
                plot(angleSpan,cosFit,'Color',pcmap(posture,:),'LineWidth',1.5)
            else
                plot(angleSpan,cosFit,'--','Color',pcmap(posture,:),'LineWidth',1.5)
            end
            hold on;
            plot(targetAngles,avgFR,'.','Color',pcmap(posture,:),'MarkerSize',5);
            xticks([]); yticks([])        
            if posture == postureList(end)
                yl = ylim;
                xl = xlim;
                scale = 0.7;
                ylRange = yl(2)-yl(1);
                ylMid = (yl(2)+yl(1))/2;
                %lowerLabel = round(ylMid-(ylRange/2)*scale);
                upperLabel = ylMid+(ylRange/2)*scale;
                text(xl(1),upperLabel,num2str(channel),'FontSize',8)
%                 text(0,lowerLabel,num2str(lowerLabel),'FontSize',8)
%                 text(0,upperLabel,num2str(upperLabel),'FontSize',8)
%                 text(xl(1)+(xl(2)-xl(1))*.7,upperLabel,['Ch' ,num2str(channel)],'FontSize',8)
            end   
            channelInd = channelInd + 1;
        end
        postureInd = postureInd + 1;
        postureListInd = postureListInd + 1;     
    end
    
    
    %% Plot tuning curves for select channels
    numPlotCh = numel(plotChannelList);
    sigThreshold = 0.10;

    
    targetAngles = (targetList-1)*45;

    fs = 14;
    
    for channel = plotChannelList
       f = figure; f.Position = [100 100 250 175]; hold on;
       for posture = postureList
            %Get tuning data
            tuningData = postureTuningData([postureTuningData.posture]==posture).tuningData; 
            avgFR = tuningData(channel).means;
            stdFR = tuningData(channel).SD;
            PD = tuningData(channel).PD;
            MD = tuningData(channel).MD;
            b0 = tuningData(channel).b0;
            p = tuningData(channel).p;
            Bfit = [b0;MD];
            angleSpan = [0:45:315]';
            x = [ones(8,1),cosd(angleSpan-PD)];
            cosFit = x*Bfit;
            %Plot
            if p < sigThreshold
                plot(angleSpan,cosFit,'Color',pcmap(posture,:),'LineWidth',2)
            else
                plot(angleSpan,cosFit,'--','Color',pcmap(posture,:),'LineWidth',1.5)
            end
            %plot(targetAngles,avgFR,'.','Color',pcmap(posture,:),'MarkerSize',5);

       end 
        %Adjust axes/labels
        yl = ylim;
        xl = xlim;
        ax = gca;
        ax.XLim = [0,315];
        xticks([0 315])
        xlabel(['Target Angle (','\circ',')'])
        set(gca, 'TickDir', 'out')
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        %Save
        if saveFig
            saveas(gcf,fullfile(saveDir,['Ch',num2str(channel),'_Reaching_TC.svg']));
        end
    end
    
    %Full Color spectrum
posture2color = [1,3,5];
    
        for channel = plotChannelList
       f = figure; f.Position = [100 100 250 175]; hold on;
       for posture = postureList
            %Get tuning data
            tuningData = postureTuningData([postureTuningData.posture]==posture).tuningData; 
            avgFR = tuningData(channel).means;
            stdFR = tuningData(channel).SD;
            PD = tuningData(channel).PD;
            MD = tuningData(channel).MD;
            b0 = tuningData(channel).b0;
            p = tuningData(channel).p;
            Bfit = [b0;MD];
            angleSpan = [0:45:315]';
            x = [ones(8,1),cosd(angleSpan-PD)];
            cosFit = x*Bfit;
            %Plot
            if p < sigThreshold
                plot(angleSpan,cosFit,'Color',pcmap(posture2color(posture),:),'LineWidth',2)
            else
                plot(angleSpan,cosFit,'--','Color',pcmap(posture2color(posture),:),'LineWidth',1.5)
            end
            %plot(targetAngles,avgFR,'.','Color',pcmap(posture,:),'MarkerSize',5);

       end 
        %Adjust axes/labels
        yl = ylim;
        xl = xlim;
        ax = gca;
        ax.XLim = [0,315];
        xticks([0 315])
        xlabel(['Target Angle (','\circ',')'])
        set(gca, 'TickDir', 'out')
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        %Save
        if saveFig
            saveas(gcf,fullfile(saveDir,['Ch',num2str(channel),'_Reaching_TC_fullSpec.svg']));
        end
    end
    %% Plot PSTH for select channels
    numPlotCh = numel(plotChannelList);
    sigThreshold = 0.10;
    
    targetAngles = (targetList-1)*45;

    fs = 14;
    target = 6;
    for channel = plotChannelList
       f = figure; f.Position = [100 100 250 175]; hold on;
       for target = [1,5]
           for posture = postureList
                %Get PSTH
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(:,channel);
                timestamps = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
                plot(timestamps,traj,'Color',pcmap(posture,:),'LineWidth',1.5)
           end 
       end
%         %Adjust axes/labels
%         yl = ylim;
%         xl = xlim;
%         ax = gca;
%         ax.XLim = [0,315];
%         xticks([0 315])
%         xlabel(['Target Angle (','\circ',')'])
%         set(gca, 'TickDir', 'out')
%         set(gca,'fontname','arial')
%         set(gca,'fontsize',fs)
%         %Save
%         if saveFig
%             saveas(gcf,fullfile(saveDir,['Ch',num2str(channel),'_TC.svg']));
%         end
    end