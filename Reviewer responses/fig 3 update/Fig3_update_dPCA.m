clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 3 update dPCA';
    %saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S5 - different joints';
    set(0, 'DefaultFigureRenderer', 'painters');
     
%% Set Parameters
   numPCsToKeep = 20;     %Num PCs to project data into before analysis
   numPdims = 2;
   numTdims = 2;
   numCVReps = 100;   %Num CV redraws
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
        'vPPMean',[],'vPPCVDist',[],'vPRCVDist',[],...
        'vPTMean',[],'vPTCVDist',[],...
        'vTTMean',[],'vTTCVDist',[],'vTRCVDist',[],...
        'vTPMean',[],'vTPCVDist',[],...
        'numSkippedFolds',[]);
    structInd = 1;
    
    %Run loop
    for datasetList = bciDatasetList%{'E20210901'}%reachDatasetList%{'E20200316','N20171215','R20201020'}%bciDatasetList
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
        [minNumCondTrials] = getMinNumCondTrials(trajStruct);
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        
        %% Compute quantities for each CV draw
        
        % Preallocate 
        vPPCVDist = NaN(1,numCVReps); vPTCVDist = NaN(1,numCVReps); 
        vTTCVDist = NaN(1,numCVReps); vTPCVDist = NaN(1,numCVReps);
        vPRCVDist = NaN(1,numCVReps); vTRCVDist = NaN(1,numCVReps);
        
        % Loop
        curInd = 1; %Index of projection vectors (e.g., vPP)
        numSkippedFolds = 0;
        for CVRep = 1:numCVReps
            CVRep
            % Divide into folds 
            [trajStruct1,trajStruct2] = splitDataforCV(trajStruct,minNumCondTrials,binWidth);  
            
            % Draw a random subspace
            v1 = normrnd(0,1,numPCsToKeep,1);
            v2 = normrnd(0,1,numPCsToKeep,1);
            DR = orth([v1,v2]);
            
            % Get P and T sig and dims (using dPCA)
            [pSig1,tSig1] = getPandTsig(trajStruct1,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
            [WGroup1,V,postureDimsGroup1,targetDimsGroup1,explVar,additionalVarExpl] = getPostureDPCADimsOnPCA(trajStruct1,true,postureList,numPostures,targetList,numTargets,numPCsToKeep,minNumTimestamps);
            
                pSig1Reshape = reshape(squeeze(pSig1),[numPostures*minNumTimestamps,numPCsToKeep]);
                tSig1Reshape = reshape(squeeze(tSig1),[numTargets*minNumTimestamps,numPCsToKeep]);
                
                pDimsGroup1 = WGroup1(:,postureDimsGroup1);
                [pDimsGroup1, ~] = qr(pDimsGroup1);
                tDimsGroup1 = WGroup1(:,targetDimsGroup1);
                [tDimsGroup1, ~] = qr(tDimsGroup1);
                
            [pSig2,tSig2] = getPandTsig(trajStruct2,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
            [WGroup2,V,postureDimsGroup2,targetDimsGroup2,explVar,additionalVarExpl] = getPostureDPCADimsOnPCA(trajStruct2,true,postureList,numPostures,targetList,numTargets,numPCsToKeep,minNumTimestamps);
                pSig2Reshape = reshape(squeeze(pSig2),[numPostures*minNumTimestamps,numPCsToKeep]);
                tSig2Reshape = reshape(squeeze(tSig2),[numTargets*minNumTimestamps,numPCsToKeep]);
                
                pDimsGroup2 = WGroup2(:,postureDimsGroup2);
                [pDimsGroup2, ~] = qr(pDimsGroup2);
                tDimsGroup2 = WGroup2(:,targetDimsGroup2);
                [tDimsGroup2, ~] = qr(tDimsGroup2);
            
            if (length(postureDimsGroup1) >= numPdims) && (length(postureDimsGroup2) >= numPdims)
                % Cross-projection
                vPP = NaN(1,2); vPT = NaN(1,2); 
                vTT = NaN(1,2); vTP = NaN(1,2);
                vPR = NaN(1,2); vTR = NaN(1,2);
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
                    vPR(projGroup) = trace(DR'*CP*DR)/trace(CP);
                    vTR(projGroup) = trace(DR'*CT*DR)/trace(CT);
                end

                vPPCVDist(CVRep) = mean(vPP);
                vPTCVDist(CVRep) = mean(vPT);
                vTTCVDist(CVRep) = mean(vTT);
                vTPCVDist(CVRep) = mean(vTP);
                vPRCVDist(CVRep) = mean(vPR);
                vTRCVDist(CVRep) = mean(vTR);
            else
                numSkippedFolds = numSkippedFolds + 1;
            end
        end
        
        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset; 
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).vPPCVDist = vPPCVDist; 
        resultStruct(structInd).vPTCVDist = vPTCVDist; 
        resultStruct(structInd).vTTCVDist = vTTCVDist; 
        resultStruct(structInd).vTPCVDist = vTPCVDist; 
        resultStruct(structInd).vPRCVDist = vPRCVDist;
        resultStruct(structInd).vTRCVDist = vTRCVDist;
        
        resultStruct(structInd).vPPMean = mean(vPPCVDist); 
        resultStruct(structInd).vPTMean = mean(vPTCVDist); 
        resultStruct(structInd).vTPMean = mean(vTPCVDist); 
        resultStruct(structInd).vTTMean = mean(vTTCVDist); 
        
        resultStruct(structInd).numSkippedFolds = numSkippedFolds;
        
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
        'vPPCVMean',[],'vPPCI',[],...
        'vPTCVMean',[],'vPTCI',[],...
        'vPRMean',[],'vPRCI',[],...
        'vPpVal',[],'vPTRpVal',[],...
        'vTTCVMean',[],'vTTCI',[],...
        'vTPCVMean',[],'vTPCI',[],...
        'vTRMean',[],'vTRCI',[],...
        'vTpVal',[],'vTPRpVal',[]);
    
    monkeyResultStructInd = 1;
    for monkey = monkeyList
        %Get monkey data
        tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));        
        %Combine across datasets to get means and distributions
        vPPCVDist = [tempResultStruct.vPPCVDist]*100; vPPCVDist = vPPCVDist(~isnan(vPPCVDist));
        
        vPPCVMean = mean(vPPCVDist);
        vPTCVDist = [tempResultStruct.vPTCVDist]*100; vPTCVDist = vPTCVDist(~isnan(vPTCVDist));
        
        vPTCVMean = mean(vPTCVDist);
        vPRCVDist = [tempResultStruct.vPRCVDist]*100; vPRCVDist = vPRCVDist(~isnan(vPRCVDist));
        
        vTTCVDist = [tempResultStruct.vTTCVDist]*100; vTTCVDist = vTTCVDist(~isnan(vTTCVDist));
        vTTCVMean = mean(vTTCVDist);
        vTPCVDist = [tempResultStruct.vTPCVDist]*100; vTPCVDist = vTPCVDist(~isnan(vTPCVDist));
        vTPCVMean = mean(vTPCVDist);
        vTRCVDist = [tempResultStruct.vTRCVDist]*100; vTRCVDist = vTRCVDist(~isnan(vTRCVDist));
        
        %Test if vPP significatly different from vPT (paired, two-tailed)
        vPPairedDiffDist = vPPCVDist-vPTCVDist;
        pd = fitdist(vPPairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vPpVal = 2*min([tail_1,tail_2]);
        
        %Test if vPT significantly different from random (paired, two-tailed)
        vPRPairedDiffDist = vPTCVDist-vPRCVDist;
        pd = fitdist(vPRPairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vPTRpVal = 2*min([tail_1,tail_2]);

        %pd = fitdist(vPRCVDist','Normal');
        %tail_1 = cdf('Normal',vPTCVMean,pd.mu,pd.sigma);
        %tail_2 = 1- cdf('Normal',vPTCVMean,pd.mu,pd.sigma);
        %vPTRpVal = 2*min([tail_1,tail_2]);
       
        %Test if vTT significantly different from vTP (paired, two-tailed)
        vTPairedDiffDist = vTTCVDist-vTPCVDist;
        pd = fitdist(vTPairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vTpVal = 2*min([tail_1,tail_2]);
        
        %Test if vTP significantly different from random (unpaired, two-tailed)
        vTRPairedDiffDist = vTPCVDist-vTRCVDist;
        pd = fitdist(vTRPairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vTPRpVal = 2*min([tail_1,tail_2]);        
        
        %pd = fitdist(vTRCVDist','Normal');
        %tail_1 = cdf('Normal',vTPCVMean,pd.mu,pd.sigma);
        %tail_2 = 1- cdf('Normal',vTPCVMean,pd.mu,pd.sigma);
        %vTPRpVal = 2*min([tail_1,tail_2]);
        
        %Add to monkeyResultStruct
        monkeyResultStruct(monkeyResultStructInd).monkey = monkey{1,1};
        
        monkeyResultStruct(monkeyResultStructInd).vPPCVMean = vPPCVMean;
        monkeyResultStruct(monkeyResultStructInd).vPPCI = [vPPCVMean-prctile(vPPCVDist,alpha*100/2) prctile(vPPCVDist,100-(alpha*100/2))-vPPCVMean];        
        monkeyResultStruct(monkeyResultStructInd).vPTCVMean = vPTCVMean;
        monkeyResultStruct(monkeyResultStructInd).vPTCI = [vPTCVMean-prctile(vPTCVDist,alpha*100/2) prctile(vPTCVDist,100-(alpha*100/2))-vPTCVMean];
        monkeyResultStruct(monkeyResultStructInd).vPpVal = vPpVal;
        monkeyResultStruct(monkeyResultStructInd).vPRMean = mean(vPRCVDist);
        monkeyResultStruct(monkeyResultStructInd).vPRCI = [prctile(vPRCVDist,alpha*100/2) prctile(vPRCVDist,100-(alpha*100/2))]; 
        monkeyResultStruct(monkeyResultStructInd).vPTRpVal = vPTRpVal;
        
        monkeyResultStruct(monkeyResultStructInd).vTTCVMean = vTTCVMean;       
        monkeyResultStruct(monkeyResultStructInd).vTTCI = [vTTCVMean-prctile(vTTCVDist,alpha*100/2) prctile(vTTCVDist,100-(alpha*100/2))-vTTCVMean];        
        monkeyResultStruct(monkeyResultStructInd).vTPCVMean =vTPCVMean;      
        monkeyResultStruct(monkeyResultStructInd).vTPCI = [vTPCVMean-prctile(vTPCVDist,alpha*100/2) prctile(vTPCVDist,100-(alpha*100/2))-vTPCVMean];        
        monkeyResultStruct(monkeyResultStructInd).vTpVal = vTpVal;
        monkeyResultStruct(monkeyResultStructInd).vTRMean = mean(vTRCVDist);
        monkeyResultStruct(monkeyResultStructInd).vTRCI = [prctile(vTRCVDist,alpha*100/2) prctile(vTRCVDist,100-(alpha*100/2))];       
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
       vPPCVMean = tempPlotStruct.vPPCVMean;
       vPPCI = tempPlotStruct.vPPCI;
       vPTCVMean = tempPlotStruct.vPTCVMean;
       vPTCI = tempPlotStruct.vPTCI;
       
       bar(monkeyInd,vPPCVMean,'FaceColor',grey,'EdgeColor',black);
       hold on;
       errorbar(monkeyInd,vPPCVMean,vPPCI(1),vPPCI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0);      
       bar(monkeyInd+offset,vPTCVMean,'FaceColor',white,'EdgeColor',black);
       errorbar(monkeyInd+offset,vPTCVMean,vPTCI(1),vPTCI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'Capsize',0);
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
       vTTCVMean = tempPlotStruct.vTTCVMean;
       vTTCI = tempPlotStruct.vTTCI;
       vTPCVMean = tempPlotStruct.vTPCVMean;
       vTPCI = tempPlotStruct.vTPCI;
       
       bar(monkeyInd,vTPCVMean,'FaceColor',grey,'EdgeColor',black);
       errorbar(monkeyInd,vTPCVMean,vTPCI(1),vTPCI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0);
       hold on;
       
       bar(monkeyInd+offset,vTTCVMean,'FaceColor',white,'EdgeColor',black);
       errorbar(monkeyInd+offset,vTTCVMean,vTTCI(1),vTTCI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0 );

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
    
    
%% Local function
    function [W,V,postureDims,targetDims,explVar,additionalVarExpl] = getPostureDPCADimsOnPCA(trajStruct,regularize,postureList,numPostures,targetList,numTargets,numChannels,minNumTimestamps)

        %Applies dPCA and returns the posture and target dims
        %used for visualizations in the paper. explVar is ordered: posture
        %1, posture 2, target 1, target 2

        %% Form X and permute to match dPCA paper
        [X] = getX(trajStruct,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);    
        Xdpca = permute(X,[4,3,2,1]);

        %% dPCA Parameters  
        numDPCs = numChannels;
        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        margNames = {'P', 'T', 'CI', 'PTI'};
        margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
        N = numChannels; P = numPostures; T = numTargets;
        time = trajStruct(1).avgZSmoothFR.timestamps(1:minNumTimestamps);
        timeEvents = [];

        %% Run dPCA
        %Regularalized version
        if regularize
            %Get Xfull and permute
            [Xfull] = getXFull(trajStruct,'allPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets); 
            XdpcaFull = permute(Xfull,[4,3,2,1,5]);
            %Set up 'numOfTrials' for use with dpca_optimizeLambda
            numOfTrials = zeros(size(Xdpca,1),size(Xdpca,2),size(Xdpca,3));
            for i = 1:size(Xdpca,1)
                for j = 1:size(Xdpca,2)
                    for k = 1:size(Xdpca,3)
                        numOfTrials(i,j,k) = sum(~isnan(XdpcaFull(i,j,k,1,:)));
                    end
                end
            end
            %Compute optimal lambda
            optimalLambda = dpca_optimizeLambda(Xdpca, XdpcaFull, numOfTrials, ...
                'numComps', numChannels,...
                'combinedParams', combinedParams, ...
                'simultaneous', true, ...
                'numRep', 10, ...  % increase this number to ~10 for better accuracy
                'display', "no", ...
                'filename', 'tmp_optimalLambdas.mat');
            %Compute Cnoise
            Cnoise = dpca_getNoiseCovariance(Xdpca, ...
                XdpcaFull, numOfTrials, 'simultaneous', true);
            % Compute dPCs
            [W,V,whichMarg] = dpca(Xdpca, numDPCs, ...
                'combinedParams', combinedParams, ...
                'lambda', optimalLambda, ...
                'Cnoise', Cnoise);

        %Unregularized version
        else
            [W,V,whichMarg] = dpca(Xdpca, numDPCs, ...
                'combinedParams', combinedParams);
        end

        %% Get variance explained by dPCs
        explVar = dpca_explainedVariance(Xdpca, W, V, ...
            'combinedParams', combinedParams);
        cumulativeDPCA = explVar.cumulativeDPCA;
        additionalVarExpl = [cumulativeDPCA(1) diff(cumulativeDPCA)];

        %% Get IDs of dPCs pertaining to posture, target, and CI (unused)
        postureMargID = find(strcmp(margNames,'P'));
        targetMargID = find(strcmp(margNames,'T'));
        CIMargID = find(strcmp(margNames,'CI'));

        postureDims = find(whichMarg==postureMargID,2);
        targetDims = find(whichMarg==targetMargID,2);
        CIDims = find(whichMarg==CIMargID,2);


    end


