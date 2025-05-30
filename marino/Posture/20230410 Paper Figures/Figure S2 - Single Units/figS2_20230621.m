clear; clc; clf; close all

%%% CHANGES TO MAKE PER STEVE MEETING %%%%%%%%%%%%%%%%%%%%
%%% - MAKE SURE CHANGE IN PD IS SIGNED
%%% - FOR EACH COMPARISON, BOOTSTRAP ESTIMATE OF TUNING CURVES, GET CHANGE IN PD
%%% - THAT GIVES YOU 10000 (SIGNED) CHANGE IN PD ESTIMATES
%%% - ASSESS WHETHER THAT DISTRIBUTION IS SIGNIFICANTLY DIFFERENT FROM ZERO
%%%   (USING T STATISTIC)



%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\S2 - Single Units';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormaps    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = orli;
    pcmap = vertcat(orli,customRainbow(4:5,:));
    scmap = customRainbow;
    anovaCmap = customRainbow([3,4,5,6],:);
    
%% Set parameters
    sigThreshold = 0.10;
    fs = 14;
    
%% Get trajStruct
    %TCDatasetList = {'E20200316','N20171215','R20201020','E20210706','N20190226','R20200221'};
    TCDatasetList = {'N20190226'};
    
    for datasetList = TCDatasetList
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'smoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318','E20200319'}
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
                task = 'bci';
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
                task = 'bci';
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
                task = 'bci';
            %Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',100}};
                task = 'reach';
            case {'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',100}};
                task = 'reach';
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',100}};
                task = 'reach';
        end 
        
       trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
      
        % Get number of trials for each condition, numTimestamps in each trial
        numCondTraj = [];
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allSmoothFR,2);
           numCondTraj = [numCondTraj,numTraj];
        end        
        figure
            histogram(numCondTraj)
            xlabel('Number of trials')
            ylabel('Number of conditions')

        %Keep only postures with all targets
        postureList = unique([trajStruct.posture]);
        targetList = unique([trajStruct.target]); 
        keepPosture = [];
        for posture = postureList
            tempTrajStruct = trajStruct([trajStruct.posture]==posture);
            postureTargetList = [tempTrajStruct.target];
            if isequal(postureTargetList,targetList)
                keepPosture = [posture,keepPosture];
            end
        end
        trajStruct = trajStruct(ismember([trajStruct.posture],keepPosture));
        
        %Get 1 observation per trial
        for i = 1:size(trajStruct,2)
            %Individual trials
            for j = 1:size(trajStruct(i).allSmoothFR,2)
               trajStruct(i).allSmoothFR(j).trialAvg = mean(trajStruct(i).allSmoothFR(j).traj); 
            end
        end
               
        
        %Get minimum number of trials and timestamps
        numCondTraj = [];
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allSmoothFR,2);
           numCondTraj = [numCondTraj,numTraj];
        end
        [minNumTrials,i] = min(numCondTraj);

        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);

%% ANOVA and Pie Plot
        %Format Data for Anova
        dataTable = NaN(minNumTrials*numPostures,numTargets,numChannels);
        postureInd = 1;
        tableInd = 1;
        for posture = postureList
            targetInd = 1;
            for target = targetList
                allSmoothFR = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allSmoothFR;
                trialAvg = vertcat(allSmoothFR.trialAvg);
                dataTable(tableInd:tableInd+minNumTrials-1,targetInd,:) = datasample(trialAvg,minNumTrials,1,'Replace',false);
                targetInd = targetInd + 1;
            end
            tableInd = tableInd + minNumTrials;
            postureInd = postureInd + 1;
        end
        
        %Perform Anova 
        resultStruct = struct('channel',[],'dataTable',[],'p',[],'tbl',[],'stats',[]);
        structInd = 1;
        for channel = 1:numChannels
            resultStruct(structInd).channel = channel;
            resultStruct(structInd).dataTable = dataTable(:,:,channel);
            [p,tbl,stats] = anova2(resultStruct(structInd).dataTable,minNumTrials,'off');
            resultStruct(structInd).p = p;
            resultStruct(structInd).tbl = tbl;
            resultStruct(structInd).stats = stats;
            structInd = structInd + 1;
        end
        
        %Plot Pie Chart
        numTarget = 0; numPosture = 0; numBoth = 0; numNeither = 0;
        for channel = 1:numChannels
            targetP = resultStruct(channel).p(1,1);
            postureP = resultStruct(channel).p(1,2);
            if postureP < 0.05 
                if targetP < 0.05
                    numBoth = numBoth + 1;
                else
                    numPosture = numPosture + 1;
                end
            end
            if postureP > 0.05 && targetP < 0.05 
                numTarget = numTarget + 1;
            end
            if postureP > 0.05 && targetP > 0.05
                numNeither = numNeither + 1;
            end
        end
        
        f = figure; f.Position = [200 200 90 90];   
        pieValues = [numTarget,numPosture,numBoth,numNeither];
        txt = {'Target Only: ';'Posture Only: ';'Mixed: ';'Neither: '}; 
        neuronCounts = {[' (',num2str(numTarget),')']; [' (',num2str(numPosture),')'];...
            [' (',num2str(numBoth),')'];[' (',num2str(numNeither),')']};        
        p = pie(pieValues(pieValues~=0));
        colormap(anovaCmap(pieValues~=0,:));       
        pText = findobj(p,'Type','text');
        percentValues = get(pText,'String'); 
        combinedtxt = strcat(txt(pieValues~=0,1),percentValues,neuronCounts(pieValues~=0,1));        
        for i = 1:size(pText,1)
            pText(i).String = combinedtxt(i);
        end        
        %title(dataset);
        %set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,dataset,'anovaPieChart.svg'));
        end

%% Get tuningData
        %Get tuning data for each posture
        postureTuningData = struct('posture',[],'tuningData',[]);
       
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
                if  any([tempTrajStruct.target]==target)
                    targetData = tempTrajStruct([tempTrajStruct.target]==target);
                    for trial = 1:numel(targetData.allSmoothFR)
                        trialFR = mean(targetData.allSmoothFR(trial).traj);
                        for channel = 1:numCh
                            tuningData(channel).allData(trial,targetInd) = trialFR(channel);
                        end
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
        
%% Get posture PD dists 
    posturePDDist = struct('channel',[],'posture',[],'PDDist',[],'PD',[],'sigTuned',[]);
    structInd = 1;
    for channel = 1:numCh
       for posture = postureList
           %Get repeated samples
           targetInd = 1;
           samples = NaN(minNumTrials,numTargets);
           for target = targetList
                allSmoothFR = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allSmoothFR;
                trialAvg = vertcat(allSmoothFR.trialAvg);
                samples(:,targetInd) =  datasample(trialAvg(:,channel),minNumTrials,1,'Replace',false);
                targetInd = targetInd + 1;
           end
           %For each sample, get PDDist
           PDDist = NaN(1,minNumTrials);
           for i = 1:minNumTrials
                y = nan(numTargets,1); x = nan(numTargets,3);
                for targetInd = 1:numTargets
                    y(targetInd,1) = samples(i,targetInd);
                    x(targetInd,:) = [1,sind(targetAngles(targetInd)),cosd(targetAngles(targetInd))];
                end
                [B,bint,r,rint,stats] = regress(y,x);
                b0 = B(1,1); b1 = B(2,1); b2 = B(3,1); p = stats(3);
                MD = sqrt(b1.^2 + b2.^2);
                PD = atan2d(b1,b2);
                if PD < 0
                    PD = 360 - abs(PD);
                end               
                PDDist(i) = PD;
           end
           %PD is mean of PDDist
           PD = mean(PDDist);
           %sigTuned if p<sigThreshold
           tempData = postureTuningData([postureTuningData.posture]==posture).tuningData;
           p = tempData([tempData.channel]==channel).p;
           if p < sigThreshold
               sigTuned = true;
           else
               sigTuned = false;
           end           
           posturePDDist(structInd).channel = channel;
           posturePDDist(structInd).posture = posture;
           posturePDDist(structInd).PDDist = PDDist;
           posturePDDist(structInd).PD = PD;
           posturePDDist(structInd).sigTuned = sigTuned;
           structInd = structInd + 1;
       end
    end
 
%% Get across-posture changes (use PD from all data here)
        acrossPostureStruct = struct('posture',[],'delPD',[],'delMD',[],'absDelPD',[],'uAbsDelPD',[],'chList',[],'sigPDChange',[]);
        structInd = 1;
        for posture = postureList(2:end)
            %Get chList to test (sig tuned in both postures)
            chList = [];
            for channel = 1:numCh
               sigTunedP1 = posturePDDist([posturePDDist.channel]==channel & [posturePDDist.posture]==1).sigTuned;
               sigTunedCurP = posturePDDist([posturePDDist.channel]==channel & [posturePDDist.posture]==posture).sigTuned;
               if sigTunedP1 && sigTunedCurP
                   chList = [chList,channel];
               end
            end
            
            %Check for significant change in tuning. If present, mark
            %sigPDChange
            delPD = NaN(1,length(chList));
            sigPDChange = false(1,length(chList));
            for i = 1:length(chList)
               channel = chList(i);
               posture1PDDist = posturePDDist([posturePDDist.channel]==channel & [posturePDDist.posture]==1).PDDist;
               %posture1PD = posturePDDist([posturePDDist.channel]==channel & [posturePDDist.posture]==1).PD;
               curPosturePDDist = posturePDDist([posturePDDist.channel]==channel & [posturePDDist.posture]==posture).PDDist;
               %curPosturePD = posturePDDist([posturePDDist.channel]==channel & [posturePDDist.posture]==posture).PD;
               [h,p,ci,stats] = ttest2(posture1PDDist,curPosturePDDist,'Alpha',sigThreshold,'Vartype','unequal');
              
               
               tuningData = postureTuningData([postureTuningData.posture]==1).tuningData;
               posture1PD = tuningData(channel).PD;
               
               tuningData = postureTuningData([postureTuningData.posture]==posture).tuningData;
               curPosturePD = tuningData(channel).PD;
               
               
               
               delPD(i) = signedAngleDiff(posture1PD,curPosturePD);
               if h == 1
                    sigPDChange(i) = true;
               end
%                figure; hold on;
%                histogram(posture1PDDist,[0:45:360]);
%                histogram(curPosturePDDist,[0:45:360]);
%                close
            end

            acrossPostureStruct(structInd).delPD = delPD;
            acrossPostureStruct(structInd).absDelPD = abs(delPD);
            acrossPostureStruct(structInd).sigPDChange = sigPDChange;
            acrossPostureStruct(structInd).posture = posture;
            acrossPostureStruct(structInd).chList = chList;
            structInd = structInd + 1;
        end
        
         
%% Make violin plots         
  %Plot Change PD violins
        violinMat = [];
        absViolinMat = [];
        i = 1;
        for posture = postureList(2:end)
            delPD = acrossPostureStruct([acrossPostureStruct.posture]==posture).delPD;
                violinMat(i:i+length(delPD)-1,1) = posture*ones(length(delPD),1);
                violinMat(i:i+length(delPD)-1,2) = delPD';
            absDelPD = acrossPostureStruct([acrossPostureStruct.posture]==posture).absDelPD;
                absViolinMat(i:i+length(absDelPD)-1,1) = posture*ones(length(absDelPD),1);
                absViolinMat(i:i+length(absDelPD)-1,2) = absDelPD';
            i = i + length(delPD);
        end
        
        plotPostures = unique(violinMat(:,1));
        fs = 20;
        
        f= figure; f.Position = [200 200 100 100]; hold on;
        violinplot(violinMat(:,2),violinMat(:,1),'ViolinColor',pcmap(plotPostures,:),'ShowData',false);
        scatterScale = 0.3;
        postureInd = 1;
        for posture = postureList(2:end)
           delPD = acrossPostureStruct([acrossPostureStruct.posture]==posture).delPD;
           sigPDChange = acrossPostureStruct([acrossPostureStruct.posture]==posture).sigPDChange;
           sigDelPD = delPD(sigPDChange);
           insigDelPD = delPD(~sigPDChange);
           plot(postureInd+scatterScale*rand(1,length(insigDelPD))-scatterScale/2,insigDelPD,'o','MarkerSize',1,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',pcmap(posture,:));
           plot(postureInd+scatterScale*rand(1,length(sigDelPD))-scatterScale/2,sigDelPD,'.','MarkerSize',6,'Color',pcmap(posture,:));
           postureInd = postureInd + 1;
        end
        %xlabel('Posture')
        %ylabel('\Delta PD rel. to Posture 1 (deg)')
        xlim([0 7])
        %set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        set(gca,'TickDir','out')
        if saveFig
            saveas(gcf,fullfile(saveDir,dataset,'ChangePD.svg'));
        end
        
       
        
        %Violin w channels labeled
        figure; hold on;
        violinplot(violinMat(:,2),violinMat(:,1),'ViolinColor',pcmap(plotPostures,:),'ShowData',false);
        scatterScale = 0.3;
        postureInd = 1;
        for posture = postureList(2:end)
           delPD = acrossPostureStruct([acrossPostureStruct.posture]==posture).delPD;
           sigPDChange = acrossPostureStruct([acrossPostureStruct.posture]==posture).sigPDChange;
           sigDelPD = delPD(sigPDChange);
           insigDelPD = delPD(~sigPDChange);
           plot(postureInd+scatterScale*rand(1,length(insigDelPD))-scatterScale/2,insigDelPD,'o','MarkerSize',3,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',pcmap(posture,:));
           plot(postureInd+scatterScale*rand(1,length(sigDelPD))-scatterScale/2,sigDelPD,'.','MarkerSize',13,'Color',pcmap(posture,:));
           postureInd = postureInd + 1;
        end

        biggestPosture = acrossPostureStruct(end).posture;
        
        chList = acrossPostureStruct(end).chList;
        for i = 1:length(chList)
           ch = chList(i);
           posture = biggestPosture;
           delPD = acrossPostureStruct(end).delPD(i);
           sigPDChange = acrossPostureStruct(end).sigPDChange;
           if abs(delPD) > 10 & sigPDChange(i)
               text(postureInd-1+randn*0.25,delPD,num2str(ch));
           end
        end
        
        xlabel('Posture')
        ylabel('\Delta PD rel. to Posture 1 (deg)')
        xlim([0 7])
        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
%         if saveFig
%             saveas(gcf,fullfile(saveDir,dataset,'ChangePD_Labelled.fig'));
%         end

%% Tuning curves
%         %Plot All TC's
%         maxChangePD = zeros(1,numCh);
%         numChannelGroups = ceil(numCh/64);
%         for channelGroup = 1:numChannelGroups
%             channelList = (channelGroup-1)*64+1:(channelGroup)*64;
%             channelList = channelList(channelList <= numCh);
%             f = figure; f.Position = [10 10 700 650];
%             [ha, pos] = tight_subplot(8,8,0,0,0);
%             postureInd = 1;
%             postureListInd = 1;
%             targetAngles = (targetList-1)*45;
%             for posture = postureList
%                 %Get tuning data
%                 tuningData = postureTuningData([postureTuningData.posture]==posture).tuningData;
%                 %Plot TC's
%                 for channel = channelList
%                     row = floor(channel-(channelList(1)-1)/8)+1;
%                     col = channel-(channelList(1)-1) - 8*(row-1); 
%                     axes(ha((row-1)*8 + col));
%                     avgFR = tuningData(channel).means;
%                     stdFR = tuningData(channel).SD;
%                     PD = tuningData(channel).PD;
%                     MD = tuningData(channel).MD;
%                     b0 = tuningData(channel).b0;
%                     p = tuningData(channel).p;
%                     Bfit = [b0;MD];
%                     angleSpan = [0:45:315]';
%                     x = [ones(8,1),cosd(angleSpan-PD)];
%                     cosFit = x*Bfit;
%                     if p < sigThreshold
%                         plot(angleSpan,cosFit,'Color',pcmap(posture,:),'LineWidth',1.5)
%                     else
%                         plot(angleSpan,cosFit,'--','Color',pcmap(posture,:),'LineWidth',1.5)
%                     end
%                     hold on;
%                     plot(targetAngles,avgFR,'.','Color',pcmap(posture,:),'MarkerSize',5);
%                     xticks([]); yticks([])        
%                     if posture == postureList(end)
%                         yl = ylim;
%                         xl = xlim;
%                         scale = 0.7;
%                         ylRange = yl(2)-yl(1);
%                         ylMid = (yl(2)+yl(1))/2;
%                         upperLabel = ylMid+(ylRange/2)*scale;
%                         lowerLabel = ylMid-(ylRange/2)*scale;
%                         text(xl(1),upperLabel,num2str(channel),'FontSize',8)
%                         
%                         channelDelPD = [];
%                         for i = 1:size(acrossPostureStruct,2)
%                             channelDelPD(i) = acrossPostureStruct(i).delPD(channel);
%                         end
%                         if any(~isnan(channelDelPD))
%                            [~,maxChannelDelPDInd] = max(abs(channelDelPD(~isnan(channelDelPD))));
%                            maxChannelDelPD = round(channelDelPD(maxChannelDelPDInd));
%                            text(xl(1),lowerLabel,['\DeltaPD ',num2str(maxChannelDelPD)],'FontSize',8)
%                         end
%                         
%         %                 text(0,lowerLabel,num2str(lowerLabel),'FontSize',8)
%         %                 text(0,upperLabel,num2str(upperLabel),'FontSize',8)
%         %                 text(xl(1)+(xl(2)-xl(1))*.7,upperLabel,['Ch' ,num2str(channel)],'FontSize',8)
%                     end         
%                 end
%                 postureInd = postureInd + 1;
%                 postureListInd = postureListInd + 1;     
%             end     
%             set(gca,'fontname','arial');
%             if saveFig
%                 saveas(gcf,fullfile(saveDir,'All TC',[dataset,'TC_ChGroup',num2str(channelGroup),'.svg']));
%             end
%         end

        %Plot select TC's
         switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318','E20200319'}
                chList = [17,11,80];
            case {'N20171215','N20180221'}
                chList = [23,56];
             case {'R20201020'}
                 chList = [20,39];
             case {'E20210706'}
                 chList = [17,88];
             case {'N20190226'}
                 chList = [14,45,55];
             case {'R20200221'}
                 chList = [75,66];
         end
         fs = 16;
         for channel = chList
             f = figure; f.Position = [200 200 50 50];
             for posture = [1,biggestPosture]
                %Get tuning data
                tuningData = postureTuningData([postureTuningData.posture]==posture).tuningData;
                allData = tuningData(channel).allData;
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
                %plot(targetAngles,avgFR,'.','Color',pcmap(posture,:),'MarkerSize',5);
                plot(targetAngles',allData,'.','Color',pcmap(posture,:),'MarkerSize',3);
                %xticks([]); yticks([]) 
                xticks([0:90:270]);
                if posture == postureList(end)
                    yl = ylim;
                    xl = xlim;
                    scale = 0.7;
                    ylRange = yl(2)-yl(1);
                    ylMid = (yl(2)+yl(1))/2;
                    upperLabel = ylMid+(ylRange/2)*scale;
                    lowerLabel = ylMid-(ylRange/2)*scale;
                    %text(xl(1),upperLabel,num2str(channel),'FontSize',8)
                end
             end
             %xlabel('Target Angle')
             %ylabel('FR (Hz)')
             ax= gca;
             xticklabels({})
             set(gca,'fontname','arial'); set(gca,'fontsize',6)
             set(gca,'TickDir','out')
             xlim([-15 330])
             if saveFig
                saveas(gcf,fullfile(saveDir,dataset,['ch',num2str(channel),'_TC.svg']));
             end
         end
         
         
    end