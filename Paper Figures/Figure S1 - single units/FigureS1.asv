clear; clc; clf; close all

%%% CHANGES TO MAKE PER STEVE MEETING %%%%%%%%%%%%%%%%%%%%
%%% - MAKE SURE CHANGE IN PD IS SIGNED
%%% - FOR EACH COMPARISON, BOOTSTRAP ESTIMATE OF TUNING CURVES, GET CHANGE IN PD
%%% - THAT GIVES YOU 10000 (SIGNED) CHANGE IN PD ESTIMATES
%%% - ASSESS WHETHER THAT DISTRIBUTION IS SIGNIFICANTLY DIFFERENT FROM ZERO
%%%   (USING T STATISTIC)

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S1 - single units';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormaps    
    [pcmap,tcmap,rainbow] = getColorMaps();
    pcmap = vertcat(pcmap,rainbow(4:5,:));
    scmap = rainbow;
    anovaCmap = rainbow([3,4,5,6],:);
    
%% Set parameters
    sigThreshold = 0.10;
    alpha = 0.05; %For 2-way ANOVA
    numBootReps = 100; %Number of bootstrap resamples when assessing significance of delPD
    trajFields = {'smoothFR'}; %Don't z-score FRs
    dataType = 'smoothFR';
    fs = 14;
%% Main loop
    %TCDatasetList = {'E20200316','N20171215','R20201020','E20210706','N20190226','R20200221'};
    TCDatasetList = {'E20200316'};
    
    for datasetList = TCDatasetList        
        %% Get trajStruct
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);           
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset,'trajFields',trajFields);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',true);      
        %Keep only postures with all targets
        [postureList,~,targetList,~,~,~] = getTrajStructDimensions(trajStruct,'dataType',dataType);
        [trajStruct] = keepOnlyPosturesWithAllTargets(trajStruct,postureList,targetList);
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct,'dataType',dataType);      
        %Get number of trials 
        [numTrials] = getTotalNumTrials(trajStruct,'dataType',dataType);
        
        %% ANOVA and Pie Plot
        %Format data for anova
        gT = zeros(1,numTrials); %grouping variables
        gP = zeros(1,numTrials);
        FR = zeros(numTrials,numChannels);
        trialInd = 1;
        for i = 1:size(trajStruct,2)
           target = trajStruct(i).target;
           posture = trajStruct(i).posture;
           numCondTrials = size(trajStruct(i).allSmoothFR,2);
           gT(trialInd:trialInd+numCondTrials-1) = ones(1,numCondTrials)*target;
           gP(trialInd:trialInd+numCondTrials-1) = ones(1,numCondTrials)*posture;
           FR(trialInd:trialInd+numCondTrials-1,:) = vertcat(trajStruct(i).allSmoothFR.trialAvg);
           trialInd = trialInd + numCondTrials;
        end
        group = {gT,gP}; %Grouping variable 
        
        %Perform anova
        resultStruct = struct('channel',[],'p',[],'tbl',[],'stats',[],'tuning',[]);
        structInd = 1;
        for channel = 1:numChannels
            resultStruct(structInd).channel = channel;
            [p,tbl,stats] = anovan(FR(:,channel)',group,'alpha',alpha,'display','off');
            resultStruct(structInd).p = p;
            resultStruct(structInd).tbl = tbl;
            resultStruct(structInd).stats = stats;
            structInd = structInd + 1;
        end
                
        %Plot Pie Chart
        %'tuning' variable: 0 = neither; 1 = target only; 2 = posture only; 3 = both;
        numTarget = 0; numPosture = 0; numBoth = 0; numNeither = 0;
        for channel = 1:numChannels
            targetP = resultStruct(channel).p(1,1);
            postureP = resultStruct(channel).p(2,1);
            if postureP < alpha
                if targetP < alpha
                    numBoth = numBoth + 1;
                    resultStruct(channel).tuning = 3;
                else
                    numPosture = numPosture + 1;
                    resultStruct(channel).tuning = 2;
                end
            end
            if postureP > alpha && targetP < alpha
                numTarget = numTarget + 1;
                resultStruct(channel).tuning = 1;
            end
            if postureP > alpha && targetP > alpha
                numNeither = numNeither + 1;
                resultStruct(channel).tuning = 0;
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
        if saveFig
            saveas(gcf,fullfile(saveDir,dataset,'anovaPieChart.svg'));
        end

        %% Get tuningData for each unit
        %Get tuning data for each posture
        postureTuningData = struct('posture',[],'tuningData',[]);
        structInd = 1;
        for posture = postureList
            % Preallocate tuningData
            tuningData = struct('channel',nan,'allData',nan(numTrials,numTargets),'means',nan(1,numTargets),'SD',nan(1,numTargets),'MD',nan,'PD',nan,'b0',nan,'p',nan);
            tuningData = repmat(tuningData,numChannels,1);
            for i = 1:numChannels
                tuningData(i).channel = i;
            end
            % Get tuningData
            tempTrajStruct = trajStruct([trajStruct.posture]==posture);
            targetInd = 1;
            for target = targetList
                if  any([tempTrajStruct.target]==target)
                    targetData = tempTrajStruct([tempTrajStruct.target]==target);
                    for trial = 1:numel(targetData.allSmoothFR)
                        trialFR = targetData.allSmoothFR(trial).trialAvg;
                        for channel = 1:numChannels
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
            for channel = 1:numChannels
                y = nan(numTargets,1); x = nan(numTargets,3);
                for targetInd = 1:numTargets
                    y(targetInd,1) = tuningData(channel).means(1,targetInd);
                    x(targetInd,:) = [1,sind(targetAngles(targetInd)),cosd(targetAngles(targetInd))];
                end
                [B,bint,r,rint,stats] = regress(y,x);
                b0 = B(1,1); b1 = B(2,1); b2 = B(3,1); p = stats(3);
                tuningData(channel).MD = sqrt(b1.^2 + b2.^2);
                PD = atan2d(b1,b2);            
                if PD < 0 %Map [-180,180] to [0,360]
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
        
        %% Assess changes in PD and significance 
        acrossPostureStruct = struct('posture',[],'delPD',[],'absDelPD',[],'chList',[],'sigPDChange',[]);
        structInd = 1;
        for posture = postureList(2:end)
            %Get list of ch that is sigTuned in ref & comparison postures
            chList = [];
            for channel = 1:numChannels
                sigTunedRefPosture = posturePDDist([posturePDDist.channel]==channel & [posturePDDist.posture]==1).sigTuned;
                sigTunedComparisonPosture = 
                if sigTunedRefPosture && sigTunedComparisonPosture
                     chList = [chList,channel];
                end
            end 
            
            %For channels that are sigTuned in both postures, assess PD
            %change and its significance 
            delPD = NaN(1,length(chList));
            sigPDChange = false(1,length(chList));

            
            for i = 1:length(chList)
               %Get PD change 
               delPD(i) = signedAngleDiff(posture1PD,curPosturePD);
               %Get empirical dist for ref posture
            
               %Get empirical dist for comp posture
               
               %Get bootstrapped delPD estimates     
               bootDelPD = nan(1,numBootReps);
               for bootRep = 1:numBootReps
                   %Get bootstrap resample for each posture
                   
                   %Get PD estimate for each posture 
                   refBootPD = 
                   compBootPD = 
                   %Measure delPD
                   bootDelPD(bootRep) = signedAngleDiff(refBootPD,compBootPD);
               end
               
               %Test if significatly different from zero
                
                
                
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
               
               
               
               
               if h == 1
                    sigPDChange(i) = true;
               end
            end
            
            %Fill struct
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