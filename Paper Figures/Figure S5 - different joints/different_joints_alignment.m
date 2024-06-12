clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S - different joints orthogonal';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Set parameters
   numPCsToKeep = 10;     %Num PCs to project data into before analysis
   numP1dims = 1;
   numP2dims = 1;
   numBootReps = 10000;   %Num bootstrap resamples
   numRandReps = 10000; %Num random subspace draws 
   cutoffNumTraj = 6; %Num trials that must be present in a condition to keep it for analysis 
   alpha = 0.05; %Significance level
   
%% Main loop
    %Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],...
        'vP1P1Mean',[],'vP1P1BootDist',[],'vP1RDist',[],...
        'vP1P2Mean',[],'vP1P2BootDist',[],...
        'vP2P2Mean',[],'vP2P2BootDist',[],'vP2RDist',[],...
        'vP2P1Mean',[],'vP2P1BootDist',[]);
    structInd = 1;
    
    %Draw random subspaces (1D)
    randSpaces = NaN(numPCsToKeep,numP1dims,numRandReps);
    for j = 1:numRandReps
        v1 = normrnd(0,1,numPCsToKeep,1);
        DR = orth(v1);
        randSpaces(:,1,j) = DR;
    end   
    
    %Run loop
    for datasetList = {'E20210901'}
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
        [trajStruct1,trajStruct2] = splitDataforCVMaxTrials(trajStruct,binWidth);  
        
        %% Fit P1 and P2 spaces for each group, compute projection variances
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        trajStruct1_P1 = trajStruct1([trajStruct1.posture]==1 | [trajStruct1.posture]==3);
        trajStruct1_P2 = trajStruct1([trajStruct1.posture]==4 | [trajStruct1.posture]==5);
        trajStruct2_P1 = trajStruct2([trajStruct2.posture]==1 | [trajStruct2.posture]==3);
        trajStruct2_P2 = trajStruct2([trajStruct2.posture]==4 | [trajStruct2.posture]==5);
        
        postureList = [1,3]; numPostures = 2;
        [p1SigGroup1,~] = getPandTsig(trajStruct1_P1,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
            p1SigGroup1Reshape = reshape(squeeze(p1SigGroup1),[numPostures*minNumTimestamps,numPCsToKeep]);
            [p1DimsGroup1,~,~,~,explainedP1,p1SigMu] = pca(p1SigGroup1Reshape); 
        
        [p1SigGroup2,~] = getPandTsig(trajStruct2_P1,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
            p1SigGroup2Reshape = reshape(squeeze(p1SigGroup2),[numPostures*minNumTimestamps,numPCsToKeep]);
            [p1DimsGroup2,~,~,~,~,~] = pca(p1SigGroup2Reshape);   
            
        postureList = [4,5]; numPostures = 2;
        [p2SigGroup1,~] = getPandTsig(trajStruct1_P2,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
            p2SigGroup1Reshape = reshape(squeeze(p2SigGroup1),[numPostures*minNumTimestamps,numPCsToKeep]);
            [p2DimsGroup1,~,~,~,explainedP2,p2SigMu] = pca(p2SigGroup1Reshape); 
        
        [p2SigGroup2,~] = getPandTsig(trajStruct2_P2,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
            p2SigGroup2Reshape = reshape(squeeze(p2SigGroup2),[numPostures*minNumTimestamps,numPCsToKeep]);
            [p2DimsGroup2,~,~,~,~,~] = pca(p2SigGroup2Reshape); 
            

        vP1P1 = NaN(1,2); vP1P2 = NaN(1,2); vP2P1 = NaN(1,2); vP2P2 = NaN(1,2);
        for projGroup = 1:2
            if projGroup == 1
                p1SigReshape = p1SigGroup1Reshape;
                p2SigReshape = p2SigGroup1Reshape;
                p1Dims = p1DimsGroup2;
                p2Dims = p2DimsGroup2;
            elseif projGroup == 2
                p1SigReshape = p1SigGroup2Reshape;
                p2SigReshape = p2SigGroup2Reshape;
                p1Dims = p1DimsGroup1;
                p2Dims = p2DimsGroup1;
            end
            CP1 = cov(p1SigReshape);
            DP1 = p1Dims(:,1:numP1dims);
            CP2 = cov(p2SigReshape);
            DP2 = p2Dims(:,1:numP2dims);
            
            vP1P1(projGroup) = trace(DP1'*CP1*DP1)./trace(CP1);
            vP1P2(projGroup) = trace(DP2'*CP1*DP2)./trace(CP1);
            vP2P2(projGroup) = trace(DP2'*CP2*DP2)./trace(CP2);
            vP2P1(projGroup) = trace(DP1'*CP2*DP1)./trace(CP2);
        end
        resultStruct(structInd).vP1P1Mean = mean(vP1P1);
        resultStruct(structInd).vP1P2Mean = mean(vP1P2);
        resultStruct(structInd).vP2P2Mean = mean(vP2P2);
        resultStruct(structInd).vP2P1Mean = mean(vP2P1);
        
        %% Compute projection variances for random spaces
        vP1RDist = NaN(1,2*numRandReps);
        vP2RDist = NaN(1,2*numRandReps);
        curInd = 1;
        for projGroup = 1:2
            if projGroup == 1
                p1SigReshape = p1SigGroup1Reshape;
                p2SigReshape = p2SigGroup1Reshape;
            elseif projGroup == 2
                p1SigReshape = p1SigGroup2Reshape;
                p2SigReshape = p2SigGroup2Reshape;
            end
            CP1 = cov(p1SigReshape);
            CP2 = cov(p2SigReshape);
            for j = 1:numRandReps
                DR = randSpaces(:,:,j);
                vP1RDist(curInd) = trace(DR'*CP1*DR)/trace(CP1);
                vP2RDist(curInd) = trace(DR'*CP2*DR)/trace(CP2);
                curInd = curInd + 1;
            end     
        end
        resultStruct(structInd).vP1RDist = vP1RDist;
        resultStruct(structInd).vP2RDist = vP2RDist;
        
        %% For each group, bootstrap resample, get P1 & P2 sig, then project into other group's space
        % Preallocate 
        vP1P1BootDist = NaN(1,2*numBootReps); vP1P2BootDist = NaN(1,2*numBootReps); 
        vP2P2BootDist = NaN(1,2*numBootReps); vP2P1BootDist = NaN(1,2*numBootReps);
        %Loop
        curInd = 1; %Index of projection vectors (e.g., vP1P1)
        for projGroup = 1:2
            if projGroup == 1
                projStruct_P1 = trajStruct1_P1;
                projStruct_P2 = trajStruct1_P2;
                p1Dims = p1DimsGroup2;
                p2Dims = p2DimsGroup2;
            elseif projGroup == 2
                projStruct_P1 = trajStruct2_P1;
                projStruct_P2 = trajStruct2_P2;
                p1Dims = p1DimsGroup1;
                p2Dims = p2DimsGroup1;
            end
            for bootRep = 1:numBootReps
                %Subsample with replacement
                sampStruct_P1 = projStruct_P1;
                for j = 1:size(sampStruct_P1,2)
                    numTraj = size(projStruct_P1(j).allPCA,2);
                    sampInd = randsample(numTraj,numTraj,true); %True => with replacement
                    sampStruct_P1(j).allPCA = projStruct_P1(j).allPCA(sampInd);
                    [sampStruct_P1(j).avgPCA.traj,sampStruct_P1(j).avgPCA.timestamps] = getAvgTraj(sampStruct_P1(j),'PCA',binWidth);               
                end
                sampStruct_P2 = projStruct_P2;
                for j = 1:size(sampStruct_P2,2)
                    numTraj = size(projStruct_P2(j).allPCA,2);
                    sampInd = randsample(numTraj,numTraj,true); %True => with replacement
                    sampStruct_P2(j).allPCA = projStruct_P2(j).allPCA(sampInd);
                    [sampStruct_P2(j).avgPCA.traj,sampStruct_P2(j).avgPCA.timestamps] = getAvgTraj(sampStruct_P2(j),'PCA',binWidth);               
                end
                
                %Get P1 and P2 signals
                postureList = [1,3]; numPostures = 2;
                [p1Sig,~] = getPandTsig(sampStruct_P1,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
                p1SigReshape = reshape(squeeze(p1Sig),[numPostures*minNumTimestamps,numPCsToKeep]);
                               
                postureList = [4,5]; numPostures = 2;
                [p2Sig,~] = getPandTsig(sampStruct_P2,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
                p2SigReshape = reshape(squeeze(p2Sig),[numPostures*minNumTimestamps,numPCsToKeep]);             
                
                %Compute P1 and P2 projections
                CP1 = cov(p1SigReshape);
                DP1 = p1Dims(:,1:numP1dims);
                CP2 = cov(p2SigReshape);
                DP2 = p2Dims(:,1:numP2dims);
                
                vP1P1BootDist(curInd) = trace(DP1'*CP1*DP1)./trace(CP1);
                vP1P2BootDist(curInd) = trace(DP2'*CP1*DP2)./trace(CP1);
                vP2P2BootDist(curInd) = trace(DP2'*CP2*DP2)./trace(CP2);
                vP2P1BootDist(curInd) = trace(DP1'*CP2*DP1)./trace(CP2);        
                %Update curInd
                curInd = curInd + 1;
                              
            end
        end       
        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset; resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).vP1P1BootDist = vP1P1BootDist; 
        resultStruct(structInd).vP1P2BootDist = vP1P2BootDist; 
        resultStruct(structInd).vP2P2BootDist = vP2P2BootDist; 
        resultStruct(structInd).vP2P1BootDist = vP2P1BootDist; 
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
        'vP1P1Mean',[],'vP1P1CI',[],...
        'vP1P2Mean',[],'vP1P2CI',[],...
        'vP1RMean',[],'vP1RCI',[],...
        'vP1pVal',[],'vP1P2RpVal',[],...
        'vP2P2Mean',[],'vP2P2CI',[],...
        'vP2P1Mean',[],'vP2P1CI',[],...
        'vP2RMean',[],'vP2RCI',[],...
        'vP2pVal',[],'vP2P1RpVal',[]);  
    
    monkeyResultStructInd = 1;
    for monkey = monkeyList
        %Get monkey data
        tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));        
        %Combine across datasets to get means and distributions
        vP1P1Mean = mean([tempResultStruct.vP1P1Mean]*100);
        vP1P1BootDist = [tempResultStruct.vP1P1BootDist]*100;
        vP1P1BootMean = mean(vP1P1BootDist);
        vP1P2Mean = mean([tempResultStruct.vP1P2Mean]*100);
        vP1P2BootDist = [tempResultStruct.vP1P2BootDist]*100;
        vP1P2BootMean = mean(vP1P2BootDist);
        vP1RDist = [tempResultStruct.vP1RDist]*100;
        
        vP2P2Mean = mean([tempResultStruct.vP2P2Mean]*100);
        vP2P2BootDist = [tempResultStruct.vP2P2BootDist]*100;
        vP2P2BootMean = mean(vP2P2BootDist);
        vP2P1Mean = mean([tempResultStruct.vP2P1Mean]*100);
        vP2P1BootDist = [tempResultStruct.vP2P1BootDist]*100;
        vP2P1BootMean = mean(vP2P1BootDist);
        vP2RDist = [tempResultStruct.vP2RDist]*100;
        
        %Test if vP1P1 significatly different from vP1P2 (paired, two-tailed)
        vP1PairedDiffDist = vP1P1BootDist-vP1P2BootDist;
        pd = fitdist(vP1PairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vP1pVal = 2*min([tail_1,tail_2]);
        
        %Test if vP1P2 significantly different from random (unpaired, two-tailed)
        pd = fitdist(vP1RDist','Normal');
        tail_1 = cdf('Normal',vP1P2Mean,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',vP1P2Mean,pd.mu,pd.sigma);
        vP1P2RpVal = 2*min([tail_1,tail_2]);
       
        %Test if vP2P2 significantly different from vP2P1 (paired, two-tailed)
        vP2PairedDiffDist = vP2P2BootDist-vP2P1BootDist;
        pd = fitdist(vP2PairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vP2pVal = 2*min([tail_1,tail_2]);
        
        %Test if vP2P1 significantly different from random (unpaired, two-tailed)
        pd = fitdist(vP2RDist','Normal');
        tail_1 = cdf('Normal',vP2P1Mean,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',vP2P1Mean,pd.mu,pd.sigma);
        vP2P1RpVal = 2*min([tail_1,tail_2]);
        
        %Add to monkeyResultStruct
        monkeyResultStruct(monkeyResultStructInd).monkey = monkey{1,1};
        
        monkeyResultStruct(monkeyResultStructInd).vP1P1Mean = vP1P1Mean;
        monkeyResultStruct(monkeyResultStructInd).vP1P1CI = [vP1P1BootMean-prctile(vP1P1BootDist,alpha*100/2) prctile(vP1P1BootDist,100-(alpha*100/2))-vP1P1BootMean];        
        monkeyResultStruct(monkeyResultStructInd).vP1P2Mean = vP1P2Mean;
        monkeyResultStruct(monkeyResultStructInd).vP1P2CI = [vP1P2BootMean-prctile(vP1P2BootDist,alpha*100/2) prctile(vP1P2BootDist,100-(alpha*100/2))-vP1P2BootMean];
        monkeyResultStruct(monkeyResultStructInd).vP1pVal = vP1pVal;
        monkeyResultStruct(monkeyResultStructInd).vP1RMean = mean(vP1RDist);
        monkeyResultStruct(monkeyResultStructInd).vP1RCI = [prctile(vP1RDist,alpha*100/2) prctile(vP1RDist,100-(alpha*100/2))]; 
        monkeyResultStruct(monkeyResultStructInd).vP1P2RpVal = vP1P2RpVal;
        
        monkeyResultStruct(monkeyResultStructInd).vP2P2Mean = vP2P2Mean;       
        monkeyResultStruct(monkeyResultStructInd).vP2P2CI = [vP2P2BootMean-prctile(vP2P2BootDist,alpha*100/2) prctile(vP2P2BootDist,100-(alpha*100/2))-vP2P2BootMean];        
        monkeyResultStruct(monkeyResultStructInd).vP2P1Mean =vP2P1Mean;      
        monkeyResultStruct(monkeyResultStructInd).vP2P1CI = [vP2P1BootMean-prctile(vP2P1BootDist,alpha*100/2) prctile(vP2P1BootDist,100-(alpha*100/2))-vP2P1BootMean];        
        monkeyResultStruct(monkeyResultStructInd).vP2pVal = vP2pVal;
        monkeyResultStruct(monkeyResultStructInd).vP2RMean = mean(vP2RDist);
        monkeyResultStruct(monkeyResultStructInd).vP2RCI = [prctile(vP2RDist,alpha*100/2) prctile(vP2RDist,100-(alpha*100/2))];       
        monkeyResultStruct(monkeyResultStructInd).vP2P1RpVal = vP2P1RpVal;
             
        monkeyResultStructInd = monkeyResultStructInd + 1;
    end
    
    %Export monkeyResultStruct
    if saveFig
        %use last dataset to determine task
        task = 'bci';
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
    
    %Bar plot for P1 signal
    f = figure;  hold on;
    f.Position = [200 200 figWidth figHeight];
    monkeyInd = 1;
    for monkey = monkeyList
       tempPlotStruct = monkeyResultStruct(strcmp(plotMonkeyList,monkey{1,1}));
       vP1P1Mean = tempPlotStruct.vP1P1Mean;
       vP1P1CI = tempPlotStruct.vP1P1CI;
       vP1P2Mean = tempPlotStruct.vP1P2Mean;
       vP1P2CI = tempPlotStruct.vP1P2CI;
       
       bar(monkeyInd,vP1P1Mean,'FaceColor',grey,'EdgeColor',black);
       hold on;
       errorbar(monkeyInd,vP1P1Mean,vP1P1CI(1),vP1P1CI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0);      
       bar(monkeyInd+offset,vP1P2Mean,'FaceColor',white,'EdgeColor',black);
       errorbar(monkeyInd+offset,vP1P2Mean,vP1P2CI(1),vP1P2CI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'Capsize',0);
       monkeyInd = monkeyInd + 1; 
    end
    
    %Combine vP1RCI across monkeys, compute mean and CI, add to plot
    vP1RMean = mean([monkeyResultStruct.vP1RMean]);
    if numMonkeys > 1
        vP1RCI = mean(vertcat(monkeyResultStruct.vP1RCI));
    else
        vP1RCI = monkeyResultStruct(1).vP1RCI;
    end
    ax = gca;
    curXLims = ax.XLim;
    shadedErrorBar(curXLims,[vP1RMean,vP1RMean],[vP1RCI(2)-vP1RMean vP1RCI(2)-vP1RMean; vP1RMean-vP1RCI(1) vP1RMean-vP1RCI(1)],'lineprops',{'--','LineWidth',.5,'Color',[0.3 0.3 0.3]});
    set(gca,'fontname','arial'); set(gca,'fontsize',fs);
    ax.TickDir = 'out';
    yticks([0:25:100]);
    xticks([]);
    ax.YLim = [0 100];
    %ylabel('Variance Explained %');
    %xlabel('Subspace');

    if saveFig
        saveas(gcf,fullfile(saveDir,[task,'_P1SigVarExpl.svg']));
    end
    
    %Bar plot for P2 signal
    f = figure;  hold on;
    f.Position = [200 200 figWidth figHeight];
    monkeyInd = 1;
    for monkey = monkeyList
       tempPlotStruct = monkeyResultStruct(strcmp(plotMonkeyList,monkey{1,1}));
       vP2P2Mean = tempPlotStruct.vP2P2Mean;
       vP2P2CI = tempPlotStruct.vP2P2CI;
       vP2P1Mean = tempPlotStruct.vP2P1Mean;
       vP2P1CI = tempPlotStruct.vP2P1CI;
       
       bar(monkeyInd,vP2P1Mean,'FaceColor',grey,'EdgeColor',black);
       errorbar(monkeyInd,vP2P1Mean,vP2P1CI(1),vP2P1CI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0);
       hold on;
       
       bar(monkeyInd+offset,vP2P2Mean,'FaceColor',white,'EdgeColor',black);
       errorbar(monkeyInd+offset,vP2P2Mean,vP2P2CI(1),vP2P2CI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0 );

       monkeyInd = monkeyInd + 1; 
    end
    
    %Combine vP12R across monkeys, compute mean and CI, add to plot
    vP2RMean = mean([monkeyResultStruct.vP2RMean]);
    if numMonkeys > 1
        vP2RCI = mean(vertcat(monkeyResultStruct.vP2RCI));
    else
        vP2RCI = monkeyResultStruct(1).vP2RCI;
    end
    ax = gca;
    curXLims = ax.XLim;
    shadedErrorBar(curXLims,[vP2RMean,vP2RMean],[vP2RCI(2)-vP2RMean vP2RCI(2)-vP2RMean; vP2RMean-vP2RCI(1) vP2RMean-vP2RCI(1)],'lineprops',{'--','LineWidth',.5,'Color',[0.3 0.3 0.3]});
    
    set(gca,'fontname','arial'); set(gca,'fontsize',fs);
    ax.TickDir = 'out';
    yticks([0:25:100]);
    xticks([]);
    %ylabel('Variance Explained %');
    %xlabel('Subspace');
        ax.YLim = [0 100];
    if saveFig
        saveas(gcf,fullfile(saveDir,[task,'_P2SigVarExpl.svg']));
    end  
    