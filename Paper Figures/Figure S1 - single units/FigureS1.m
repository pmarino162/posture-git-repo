clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S1 - single units';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set parameters
    refPosture = 1; %Use posture 1 for all task x animals
    alpha = 0.05; %Significance level for all analyses
    numBootReps = 10000; %Number of bootstrap resamples when assessing significance of delPD
    
%% Main loop
    TCDatasetList = {'E20200316','N20171215','R20201020','E20210706','N20190226','R20200221'};
    %E20200316 1=I30 - 5=E30
    %E20210901 1=I30 3=E30 4=F30 5=eE30
    %N20171215 1=N00, 2=I45, 3=A90
    %R20201020 1=N00, 2=A90
    %R20201021 1=N00, 2=I45 
    %E20210706
    %N20190226 1=(36,-60), 2=(-66,9)
    %R20200211 1=(), 2=()
    
    %TCDatasetList = {'R20200221'};
    
    for datasetList = TCDatasetList
        %% Get trajStruct
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);           
        %Get trajStruct
        [condFields,~,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajFields = {'singleBinFR'}; %Use one big bin - no smoothing or z-scoring
        switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318'}
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            case {'N20171215','N20180221'}
                trialInclStates(1).inclStates = {{'state','Cursor Release','first',0},{'state','Target Hold','first',0}};
            case {'R20201020','R20201021'}
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
            %Reaching
            case 'E20210706'
                trialInclStates(1).inclStates = {{'state','Target Acquire','first',0},{'state','Target Hold','first',0}};
            case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
                trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Target Hold','first',0}};
            case{'R20200221','R20200222'}
                trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Target Hold','first',0}};
        end
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);      
        %Keep only postures with all targets
        [postureList,~,targetList,~,~,~] = getTrajStructDimensions(trajStruct);
        [trajStruct] = keepOnlyPosturesWithAllTargets(trajStruct,postureList,targetList);
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);      
        %Get number of trials 
        [numTrials] = getTotalNumTrials(trajStruct);
        
        %% Setup colormap (based on number of postures)
        [pcmap,tcmap,rainbow] = getColorMaps(numPostures);  
        if strcmpi(dataset,'E20210706') %For this analysis only, update Monkey E reaching color scheme. Helps with figure consistency 
            pcmap(5:7,:) = pcmap([6,7,5],:);
        end
        anovaCmap = rainbow([3,4,5,6],:);
        
        %% ANOVA and Pie Plot
        %Format data for anova
        gT = zeros(1,numTrials); %grouping variables
        gP = zeros(1,numTrials);
        FR = zeros(numTrials,numChannels);
        trialInd = 1;
        for i = 1:size(trajStruct,2)
           target = trajStruct(i).target;
           posture = trajStruct(i).posture;
           numCondTrials = size(trajStruct(i).allSingleBinFR,2);
           gT(trialInd:trialInd+numCondTrials-1) = ones(1,numCondTrials)*target;
           gP(trialInd:trialInd+numCondTrials-1) = ones(1,numCondTrials)*posture;
           FR(trialInd:trialInd+numCondTrials-1,:) = vertcat(trajStruct(i).allSingleBinFR.traj);
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
                
        %Count each type of neuron
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
        
        %Export anova resultStruct
        if saveFig
            save(fullfile(saveDir,dataset,'anovaResult.mat'),'resultStruct');
        end
        
        %Plot Pie Chart
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
            tuningData = struct('channel',nan,'allData',nan(numTrials,numTargets),'means',nan(1,numTargets),'SD',nan(1,numTargets),'MD',nan,'PD',nan,'b0',nan,'p',nan,'sigTuned',nan);
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
                    for trial = 1:numel(targetData.allSingleBinFR)
                        trialFR = targetData.allSingleBinFR(trial).traj;
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
            for channel = 1:numChannels
                targetMeans = tuningData(channel).means;
                [PD,MD,b0,p] = fitTC(targetMeans);
                tuningData(channel).PD= PD;
                tuningData(channel).MD = MD;
                tuningData(channel).b0 = b0;
                tuningData(channel).p = p;
                if p < alpha
                    tuningData(channel).sigTuned = true;
                else
                    tuningData(channel).sigTuned = false;
                end
            end
            %Fill Struct
            postureTuningData(structInd).posture = posture;
            postureTuningData(structInd).tuningData = tuningData;
            structInd = structInd + 1;
        end
        
        %% Assess changes in PD and significance 
        acrossPostureStruct = struct('compPosture',[],'delPD',[],'absDelPD',[],'chList',[],'sigPDChange',[]);
        structInd = 1;
        for compPosture = postureList(2:end)
            compPosture
            %Get list of ch that is sigTuned in ref & comparison postures
            chList = [];
            for channel = 1:numChannels
                sigTunedRefPosture = postureTuningData([postureTuningData.posture]==refPosture).tuningData(channel).sigTuned;
                sigTunedCompPosture = postureTuningData([postureTuningData.posture]==compPosture).tuningData(channel).sigTuned;
                if sigTunedRefPosture && sigTunedCompPosture
                     chList = [chList,channel];
                end
            end 
            
            %For channels that are sigTuned in both postures, assess PD
            %change and its significance 
            delPD = NaN(1,length(chList));
            sigPDChange = false(1,length(chList));          
            for i = 1:length(chList)
               channel = chList(i);
               %Get PD change 
               refPosturePD = postureTuningData([postureTuningData.posture]==refPosture).tuningData(channel).PD;
               compPosturePD = postureTuningData([postureTuningData.posture]==compPosture).tuningData(channel).PD;
               delPD(i) = signedAngleDiff(refPosturePD,compPosturePD);
               %Get empirical dist for ref posture
               refEmpDist = postureTuningData([postureTuningData.posture]==refPosture).tuningData(channel).allData;
               %Get empirical dist for comp posture
               compEmpDist = postureTuningData([postureTuningData.posture]==compPosture).tuningData(channel).allData;
               %Get bootstrapped delPD estimates     
               bootDelPD = nan(1,numBootReps);               
               for bootRep = 1:numBootReps
                   %Get bootstrap resample for each posture
                   refBootResample = nan(size(refEmpDist));
                   compBootResample = nan(size(compEmpDist));
                   for target = 1:8
                      numTargObs = sum(~isnan(refEmpDist(:,target)));
                      refBootResample(1:numTargObs,target) = datasample(refEmpDist(1:numTargObs,target),numTargObs);
                      numTargObs = sum(~isnan(compEmpDist(:,target)));
                      compBootResample(1:numTargObs,target) = datasample(compEmpDist(1:numTargObs,target),numTargObs);
                   end                  
                   %Compute target means for each posture 
                   refTargetMeans = nanmean(refBootResample);
                   compTargetMeans = nanmean(compBootResample);
                   %Get PD estimate for each posture 
                   [refBootPD,~,~,~] = fitTC(refTargetMeans); 
                   [compBootPD,~,~,~] = fitTC(compTargetMeans); 
                   %Measure delPD
                   bootDelPD(bootRep) = signedAngleDiff(refBootPD,compBootPD);
               end              
               %Test if significatly different from zero (two-tailed)
               CI = [prctile(bootDelPD,(alpha*100/2)),prctile(bootDelPD,100-(alpha*100/2))];
               if ~(CI(1)<0 && CI(2)>0) %if CI does not contain zero
                   sigPDChange(i) = true;
               end
            end
            
            %Fill struct
            acrossPostureStruct(structInd).delPD = delPD;
            acrossPostureStruct(structInd).absDelPD = abs(delPD);
            acrossPostureStruct(structInd).sigPDChange = sigPDChange;
            acrossPostureStruct(structInd).compPosture = compPosture;
            acrossPostureStruct(structInd).chList = chList;
            structInd = structInd + 1;            
        end

        %% Make violin plots         
        %Plot Change PD violins
        violinMat = []; 
        i = 1;
        for compPosture = postureList(2:end)
            delPD = acrossPostureStruct([acrossPostureStruct.compPosture]==compPosture).delPD;
                violinMat(i:i+length(delPD)-1,1) = compPosture*ones(length(delPD),1);
                violinMat(i:i+length(delPD)-1,2) = delPD';
            i = i + length(delPD);
        end       
        plotPostures = unique(violinMat(:,1));
        fs = 20;      
        
        f= figure; f.Position = [200 200 100 100]; hold on;
        violinplot(violinMat(:,2),violinMat(:,1),'ViolinColor',pcmap(plotPostures,:),'ShowData',false,'ShowBox',false,'ShowMedian',false,'ShowWhiskers',false);
        scatterScale = 0.3;
        compPostureInd = 1;
        for compPosture = postureList(2:end)
           delPD = acrossPostureStruct([acrossPostureStruct.compPosture]==compPosture).delPD;
           sigPDChange = acrossPostureStruct([acrossPostureStruct.compPosture]==compPosture).sigPDChange;
           sigDelPD = delPD(sigPDChange);
           insigDelPD = delPD(~sigPDChange);
           plot(compPostureInd+scatterScale*rand(1,length(insigDelPD))-scatterScale/2,insigDelPD,'o','MarkerSize',1,'MarkerEdgeColor',[.25 .25 .25],'LineWidth',0.001);
           plot(compPostureInd+scatterScale*rand(1,length(sigDelPD))-scatterScale/2,sigDelPD,'o','MarkerSize',2,'MarkerFaceColor',pcmap(compPosture,:),'MarkerEdgeColor',[0 0 0],'LineWidth',0.01);
           compPostureInd = compPostureInd + 1;
        end
        xlim([0 7])
        set(gca,'TickDir','out')
        if saveFig
            saveas(gcf,fullfile(saveDir,dataset,'ChangePD.svg'));
        end
        
        %Violin w channels labeled
        figure; hold on;
        violinplot(violinMat(:,2),violinMat(:,1),'ViolinColor',pcmap(plotPostures,:),'ShowData',false,'ShowBox',false,'ShowMedian',false,'ShowWhiskers',false);
        scatterScale = 0.3;
        compPostureInd = 1;
        for compPosture = postureList(2:end)
           delPD = acrossPostureStruct([acrossPostureStruct.compPosture]==compPosture).delPD;
           sigPDChange = acrossPostureStruct([acrossPostureStruct.compPosture]==compPosture).sigPDChange;
           sigDelPD = delPD(sigPDChange);
           insigDelPD = delPD(~sigPDChange);
           plot(compPostureInd+scatterScale*rand(1,length(insigDelPD))-scatterScale/2,insigDelPD,'o','MarkerSize',3,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',pcmap(compPosture,:));
           plot(compPostureInd+scatterScale*rand(1,length(sigDelPD))-scatterScale/2,sigDelPD,'.','MarkerSize',13,'Color',pcmap(compPosture,:));
           compPostureInd = compPostureInd + 1;
        end
        biggestPosture = acrossPostureStruct(end).compPosture;        
        chList = acrossPostureStruct(end).chList;
        for i = 1:length(chList)
           ch = chList(i);
           compPosture = biggestPosture;
           delPD = acrossPostureStruct(end).delPD(i);
           sigPDChange = acrossPostureStruct(end).sigPDChange;
           if abs(delPD) > 10 & sigPDChange(i)
               text(compPostureInd-1+randn*0.25,delPD,num2str(ch));
           end
        end       
        xlabel('Posture')
        ylabel('\Delta PD rel. to Posture 1 (deg)')
        xlim([0 7])
        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,dataset,'ChangePD_Labelled.fig'));
            saveas(gcf,fullfile(saveDir,dataset,'ChangePD_Labelled.svg'));
        end

        %Plot select TC's
         switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318','E20200319'}
                chList = [50,80,82,24,66,61];
            case {'N20171215','N20180221'}
                chList = [23,56,7,12,30,55,76];
             case {'R20201020'}
                 chList = [11,3,43,20,39,58,61,1];
             case {'E20210706'}
                 chList = [17,56,49,88,97];
             case {'N20190226'}
                 chList = [4,22,14,45,55,31,26,33];
             case {'R20200221'}
                 chList = [19,59,83,18,77,87,141,133,97,132,111];
         end
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
                if p < alpha
                    plot(angleSpan,cosFit,'Color',pcmap(posture,:),'LineWidth',1.5)
                else
                    plot(angleSpan,cosFit,'--','Color',pcmap(posture,:),'LineWidth',1.5)
                end
                hold on;
                plot(angleSpan',allData,'.','Color',pcmap(posture,:),'MarkerSize',3);
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