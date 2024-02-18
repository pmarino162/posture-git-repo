clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 3';
    %saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S5 - different joints';
    set(0, 'DefaultFigureRenderer', 'painters');
     
%% Set Parameters
   numPCsToKeep = 10;     %Num PCs to project data into before analysis
   numPdims = 2;
   numTdims = 2;
   numBootReps = 10000;   %Num bootstrap resamples
   numRandReps = 10000; %Num random subspace draws 
   cutoffNumTraj = 10; %Num trials that must be present in a condition to keep it for analysis 
   alpha = 0.05; %Significance level
   
%% Datasets to include in analysis 
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','E20210901','N20171215','N20180221','R20201020','R20201021'};         
    isoDatasetList = {'E20200116','E20200117','E20200120'};
    
%% Main loop   
    %Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],...
        'vPPMean',[],'vPPBootDist',[],'vPRDist',[],...
        'vPTMean',[],'vPTBootDist',[],...
        'vTTMean',[],'vTTBootDist',[],'vTRDist',[],...
        'vTPMean',[],'vTPBootDist',[]);
    structInd = 1;
    
    %Draw random subspaces
    randSpaces = NaN(numPCsToKeep,numPdims,numRandReps);
    for j = 1:numRandReps
        v1 = normrnd(0,1,numPCsToKeep,1);
        v2 = normrnd(0,1,numPCsToKeep,1);
        DR = orth([v1,v2]);
        randSpaces(:,:,j) = DR;
    end      
    
    %Run loop
    for datasetList = {'E20210901'}%reachDatasetList%{'E20200316','N20171215','R20201020'}%bciDatasetList
        %% Set up trajStruct
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);            
        %Project all data down to top PCs
        [allTraj] = collectAllAvgTraj(trajStruct);
        [trajStruct] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep);
        %Remove any conditions for which there weren't enough trials
        [numCondTrials] = getNumCondTrials(trajStruct,'showHist',false);         
        trajStruct = trajStruct(numCondTrials >= cutoffNumTraj);       
        %Keep only postures with all targets; get trajStructDims
        [postureList,~,targetList,~,~,~] = getTrajStructDimensions(trajStruct);
        [trajStruct] = keepOnlyPosturesWithAllTargets(trajStruct,postureList,targetList);
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);

 
        %% Split data into two groups for cross-validation
        [minNumCondTrials] = getMinNumCondTrials(trajStruct);
        [trajStruct1,trajStruct2] = splitDataforCV(trajStruct,minNumCondTrials,binWidth);  
        
        %% Fit P and T spaces for each group, compute projection variances
        [minNumTimestamps] = getMinNumTimestamps(trajStruct);    
        [pSig1,tSig1] = getPandTsig(trajStruct1,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
            pSig1Reshape = reshape(squeeze(pSig1),[numPostures*minNumTimestamps,numPCsToKeep]);
            tSig1Reshape = reshape(squeeze(tSig1),[numTargets*minNumTimestamps,numPCsToKeep]);
            [pDimsGroup1,~,~,~,explainedP,pSigMu] = pca(pSig1Reshape); 
            [tDimsGroup1,~,~,~,explainedT,tSigMu] = pca(tSig1Reshape); 
        [pSig2,tSig2] = getPandTsig(trajStruct2,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
            pSig2Reshape = reshape(squeeze(pSig2),[numPostures*minNumTimestamps,numPCsToKeep]);
            tSig2Reshape = reshape(squeeze(tSig2),[numTargets*minNumTimestamps,numPCsToKeep]);
            [pDimsGroup2,~,~,~,~,~] = pca(pSig2Reshape); 
            [tDimsGroup2,~,~,~,~,~] = pca(tSig2Reshape); 
        vPP = NaN(1,2); vPT = NaN(1,2); vTT = NaN(1,2); vTP = NaN(1,2);
        for projGroup = 1:2
            if projGroup == 1
                pSigReshape = pSig1Reshape;
                tSigReshape = tSig1Reshape;
                pDims = pDimsGroup2;
                tDims = tDimsGroup2;
            elseif projGroup == 2
                pSigReshape = pSig2Reshape;
                tSigReshape = tSig2Reshape;
                pDims = pDimsGroup1;
                tDims = tDimsGroup1;
            end
            CP = cov(pSigReshape);
            DP = pDims(:,1:numPdims);
            CT = cov(tSigReshape);
            DT = tDims(:,1:numTdims);
            vPP(projGroup) = trace(DP'*CP*DP)./trace(CP);
            vPT(projGroup) = trace(DT'*CP*DT)./trace(CP);
            vTT(projGroup) = trace(DT'*CT*DT)./trace(CT);
            vTP(projGroup) = trace(DP'*CT*DP)./trace(CT);
        end
        resultStruct(structInd).vPPMean = mean(vPP);
        resultStruct(structInd).vPTMean = mean(vPT);
        resultStruct(structInd).vTTMean = mean(vTT);
        resultStruct(structInd).vTPMean = mean(vTP);
        
        %% Compute projection variance for random spaces
        vPRDist = NaN(1,2*numRandReps);
        vTRDist = NaN(1,2*numRandReps);
        curInd = 1;
        for projGroup = 1:2
            if projGroup == 1
                pSigReshape = pSig1Reshape;
                tSigReshape = tSig1Reshape;
            elseif projGroup == 2
                pSigReshape = pSig2Reshape;
                tSigReshape = tSig2Reshape;
            end
            CP = cov(pSigReshape);
            CT = cov(tSigReshape);
            for j = 1:numRandReps
                DR = randSpaces(:,:,j);
                vPRDist(curInd) = trace(DR'*CP*DR)/trace(CP);
                vTRDist(curInd) = trace(DR'*CT*DR)/trace(CT);
                curInd = curInd + 1;
            end     
        end
        resultStruct(structInd).vPRDist = vPRDist;
        resultStruct(structInd).vTRDist = vTRDist;
        
        %% For each group, bootstrap resample, get P & T sig, then project into other group's space
        % Preallocate 
        vPPBootDist = NaN(1,2*numBootReps); vPTBootDist = NaN(1,2*numBootReps); 
        vTTBootDist = NaN(1,2*numBootReps); vTPBootDist = NaN(1,2*numBootReps);
        %Loop
        curInd = 1; %Index of projection vectors (e.g., vPP)
        for projGroup = 1:2
            if projGroup == 1
                projStruct = trajStruct1;
                pDims = pDimsGroup2;
                tDims = tDimsGroup2;
            elseif projGroup == 2
                projStruct = trajStruct2;
                pDims = pDimsGroup1;
                tDims = tDimsGroup1;
            end
            for bootRep = 1:numBootReps
                %Subsample with replacement
                sampStruct = projStruct;
                for j = 1:size(sampStruct,2)
                    numTraj = size(projStruct(j).allPCA,2);
                    sampInd = randsample(numTraj,numTraj,true); %True => with replacement
                    sampStruct(j).allPCA = projStruct(j).allPCA(sampInd);
                    [sampStruct(j).avgPCA.traj,sampStruct(j).avgPCA.timestamps] = getAvgTraj(sampStruct(j),'PCA',binWidth);               
                end
                %Get P and T signals
                [pSig,tSig] = getPandTsig(sampStruct,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
                pSigReshape = reshape(squeeze(pSig),[numPostures*minNumTimestamps,numPCsToKeep]);
                tSigReshape = reshape(squeeze(tSig),[numTargets*minNumTimestamps,numPCsToKeep]);
                %Compute P and T projections
                CP = cov(pSigReshape);
                DP = pDims(:,1:numPdims);
                CT = cov(tSigReshape);
                DT = tDims(:,1:numTdims);
                vPPBootDist(curInd) = trace(DP'*CP*DP)./trace(CP);
                vPTBootDist(curInd) = trace(DT'*CP*DT)./trace(CP);
                vTTBootDist(curInd) = trace(DT'*CT*DT)./trace(CT);
                vTPBootDist(curInd) = trace(DP'*CT*DP)./trace(CT);        
                %Update curInd
                curInd = curInd + 1;
            end
        end       
        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset; resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).vPPBootDist = vPPBootDist; 
        resultStruct(structInd).vPTBootDist = vPTBootDist; 
        resultStruct(structInd).vTTBootDist = vTTBootDist; 
        resultStruct(structInd).vTPBootDist = vTPBootDist; 
        structInd = structInd + 1;       
    end
    
%% Get monkey list 
    %Across datasets
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
%% For each monkey, get overall distributions, CI's, and p-values
    monkeyResultStruct = struct('monkey',[],...
        'vPPMean',[],'vPPCI',[],...
        'vPTMean',[],'vPTCI',[],...
        'vPRMean',[],'vPRCI',[],...
        'vPpVal',[],'vPTRpVal',[],...
        'vTTMean',[],'vTTCI',[],...
        'vTPMean',[],'vTPCI',[],...
        'vTRMean',[],'vTRCI',[],...
        'vTpVal',[],'vTPRpVal',[]);
    
    monkeyResultStructInd = 1;
    for monkey = monkeyList
        %Get monkey data
        tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));        
        %Combine across datasets to get means and distributions
        vPPMean = mean([tempResultStruct.vPPMean]*100);
        vPPBootDist = [tempResultStruct.vPPBootDist]*100;
        vPPBootMean = mean(vPPBootDist);
        vPTMean = mean([tempResultStruct.vPTMean]*100);
        vPTBootDist = [tempResultStruct.vPTBootDist]*100;
        vPTBootMean = mean(vPTBootDist);
        vPRDist = [tempResultStruct.vPRDist]*100;
        
        vTTMean = mean([tempResultStruct.vTTMean]*100);
        vTTBootDist = [tempResultStruct.vTTBootDist]*100;
        vTTBootMean = mean(vTTBootDist);
        vTPMean = mean([tempResultStruct.vTPMean]*100);
        vTPBootDist = [tempResultStruct.vTPBootDist]*100;
        vTPBootMean = mean(vTPBootDist);
        vTRDist = [tempResultStruct.vTRDist]*100;
        
        %Test if vPP significatly different from vPT (paired, two-tailed)
        vPPairedDiffDist = vPPBootDist-vPTBootDist;
        pd = fitdist(vPPairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vPpVal = 2*min([tail_1,tail_2]);
        
        %Test if vPT significantly different from random (unpaired, two-tailed)
        pd = fitdist(vPRDist','Normal');
        tail_1 = cdf('Normal',vPTMean,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',vPTMean,pd.mu,pd.sigma);
        vPTRpVal = 2*min([tail_1,tail_2]);
       
        %Test if vTT significantly different from vTP (paired, two-tailed)
        vTPairedDiffDist = vTTBootDist-vTPBootDist;
        pd = fitdist(vTPairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vTpVal = 2*min([tail_1,tail_2]);
        
        %Test if vTP significantly different from random (unpaired, two-tailed)
        pd = fitdist(vTRDist','Normal');
        tail_1 = cdf('Normal',vTPMean,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',vTPMean,pd.mu,pd.sigma);
        vTPRpVal = 2*min([tail_1,tail_2]);
        
        %Add to monkeyResultStruct
        monkeyResultStruct(monkeyResultStructInd).monkey = monkey{1,1};
        
        monkeyResultStruct(monkeyResultStructInd).vPPMean = vPPMean;
        monkeyResultStruct(monkeyResultStructInd).vPPCI = [vPPBootMean-prctile(vPPBootDist,alpha*100/2) prctile(vPPBootDist,100-(alpha*100/2))-vPPBootMean];        
        monkeyResultStruct(monkeyResultStructInd).vPTMean = vPTMean;
        monkeyResultStruct(monkeyResultStructInd).vPTCI = [vPTBootMean-prctile(vPTBootDist,alpha*100/2) prctile(vPTBootDist,100-(alpha*100/2))-vPTBootMean];
        monkeyResultStruct(monkeyResultStructInd).vPpVal = vPpVal;
        monkeyResultStruct(monkeyResultStructInd).vPRMean = mean(vPRDist);
        monkeyResultStruct(monkeyResultStructInd).vPRCI = [prctile(vPRDist,alpha*100/2) prctile(vPRDist,100-(alpha*100/2))]; 
        monkeyResultStruct(monkeyResultStructInd).vPTRpVal = vPTRpVal;
        
        monkeyResultStruct(monkeyResultStructInd).vTTMean = vTTMean;       
        monkeyResultStruct(monkeyResultStructInd).vTTCI = [vTTBootMean-prctile(vTTBootDist,alpha*100/2) prctile(vTTBootDist,100-(alpha*100/2))-vTTBootMean];        
        monkeyResultStruct(monkeyResultStructInd).vTPMean =vTPMean;      
        monkeyResultStruct(monkeyResultStructInd).vTPCI = [vTPBootMean-prctile(vTPBootDist,alpha*100/2) prctile(vTPBootDist,100-(alpha*100/2))-vTPBootMean];        
        monkeyResultStruct(monkeyResultStructInd).vTpVal = vTpVal;
        monkeyResultStruct(monkeyResultStructInd).vTRMean = mean(vTRDist);
        monkeyResultStruct(monkeyResultStructInd).vTRCI = [prctile(vTRDist,alpha*100/2) prctile(vTRDist,100-(alpha*100/2))];       
        monkeyResultStruct(monkeyResultStructInd).vTPRpVal = vTPRpVal;
             

        
        monkeyResultStructInd = monkeyResultStructInd + 1;
    end
    
    %Export monkeyResultStruct
    if saveFig
        %use last dataset to determine task
        if strcmpi(dataset,'R20201021')
            task = 'bci';
        elseif strcmpi(dataset,'E20200120')
            task = 'iso';
        elseif strcmpi(dataset,'R20200222')
            task = 'reach';
        end
        save(fullfile(saveDir,[task,'_monkeyResultStruct.mat']),'monkeyResultStruct');
    end  
    
%% Plot results
    white = [1 1 1];
    grey = .75.*[1 1 1];
    red = [1 0 0];
    cyan = [0 1 1];
    black = [0 0 0];    
    mcmap = vertcat(black,cyan,red);
    figHeight = 100;
    figWidth = 10;

    %Plotting parameters
    fs = 6;
    offset = 5;   
    for i = 1:numel(monkeyResultStruct)
        plotMonkeyList{i} = monkeyResultStruct(i).monkey;
    end
    
    %Bar plot for posture signal
    f = figure;  hold on;
    f.Position = [200 200 figWidth figHeight];
    monkeyInd = 1;
    for monkey = monkeyList
       tempPlotStruct = monkeyResultStruct(strcmp(plotMonkeyList,monkey{1,1}));
       vPPMean = tempPlotStruct.vPPMean;
       vPPCI = tempPlotStruct.vPPCI;
       vPTMean = tempPlotStruct.vPTMean;
       vPTCI = tempPlotStruct.vPTCI;
       
       bar(monkeyInd,vPPMean,'FaceColor',grey,'EdgeColor',black);
       hold on;
       errorbar(monkeyInd,vPPMean,vPPCI(1),vPPCI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0);      
       bar(monkeyInd+offset,vPTMean,'FaceColor',white,'EdgeColor',black);
       errorbar(monkeyInd+offset,vPTMean,vPTCI(1),vPTCI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'Capsize',0);
       monkeyInd = monkeyInd + 1; 
    end
    
    %Combine vPRCI across monkeys, compute mean and CI, add to plot
    vPRMean = mean([monkeyResultStruct.vPRMean]);
    if numMonkeys > 1
        vPRCI = mean(vertcat(monkeyResultStruct.vPRCI));
    else
        vPRCI = monkeyResultStruct(1).vPRCI;
    end
    ax = gca;
    curXLims = ax.XLim;
    shadedErrorBar(curXLims,[vPRMean,vPRMean],[vPRCI(2)-vPRMean vPRCI(2)-vPRMean; vPRMean-vPRCI(1) vPRMean-vPRCI(1)],'lineprops',{'--','LineWidth',.5,'Color',[0.3 0.3 0.3]});
    set(gca,'fontname','arial'); set(gca,'fontsize',fs);
    ax.TickDir = 'out';
    yticks([0:25:100]);
    xticks([]);
    ax.YLim = [0 100];
    %ylabel('Variance Explained %');
    %xlabel('Subspace');

    if saveFig
        saveas(gcf,fullfile(saveDir,[task,'_postureSigVarExpl.svg']));
    end
    
    %Bar plot for goal signal
    f = figure;  hold on;
    f.Position = [200 200 figWidth figHeight];
    monkeyInd = 1;
    for monkey = monkeyList
       tempPlotStruct = monkeyResultStruct(strcmp(plotMonkeyList,monkey{1,1}));
       vTTMean = tempPlotStruct.vTTMean;
       vTTCI = tempPlotStruct.vTTCI;
       vTPMean = tempPlotStruct.vTPMean;
       vTPCI = tempPlotStruct.vTPCI;
       
       bar(monkeyInd,vTPMean,'FaceColor',grey,'EdgeColor',black);
       errorbar(monkeyInd,vTPMean,vTPCI(1),vTPCI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0);
       hold on;
       
       bar(monkeyInd+offset,vTTMean,'FaceColor',white,'EdgeColor',black);
       errorbar(monkeyInd+offset,vTTMean,vTTCI(1),vTTCI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0 );

       monkeyInd = monkeyInd + 1; 
    end
    
    %Combine vTR across monkeys, compute mean and CI, add to plot
    vTRMean = mean([monkeyResultStruct.vTRMean]);
    if numMonkeys > 1
        vTRCI = mean(vertcat(monkeyResultStruct.vTRCI));
    else
        vTRCI = monkeyResultStruct(1).vTRCI;
    end
    ax = gca;
    curXLims = ax.XLim;
    shadedErrorBar(curXLims,[vTRMean,vTRMean],[vTRCI(2)-vTRMean vTRCI(2)-vTRMean; vTRMean-vTRCI(1) vTRMean-vTRCI(1)],'lineprops',{'--','LineWidth',.5,'Color',[0.3 0.3 0.3]});
    
    set(gca,'fontname','arial'); set(gca,'fontsize',fs);
    ax.TickDir = 'out';
    yticks([0:25:100]);
    xticks([]);
    %ylabel('Variance Explained %');
    %xlabel('Subspace');
        ax.YLim = [0 100];
    if saveFig
        saveas(gcf,fullfile(saveDir,[task,'_goalSigVarExpl.svg']));
    end