clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive\Documents\Posture\Paper\Reviewer responses\Analyses\Fig S5 update';
    set(0, 'DefaultFigureRenderer', 'painters');
     
%% Set Parameters
   numPCsToKeep = 10;     %Num PCs to project data into before analysis
   numP1dims = 1;
   numP2dims = 1;
   numCVReps = 10000;   %Num CV redraws
   cutoffNumTraj = 6; %Num trials that must be present in a condition to keep it for analysis 
   alpha = 0.05; %Significance level
   
    %Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],...
        'vP1P1Mean',[],'vP1P1CVDist',[],'vP1RDist',[],...
        'vP1P2Mean',[],'vP1P2CVDist',[],...
        'vP2P2Mean',[],'vP2P2CVDist',[],'vP2RDist',[],...
        'vP2P1Mean',[],'vP2P1CVDist',[]);
    structInd = 1;
    
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
        [minNumCondTrials] = getMinNumCondTrials(trajStruct);
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        
        %% Compute quantities for each CV draw
        
        % Preallocate 
        vP1P1CVDist = NaN(1,numCVReps); vP1P2CVDist = NaN(1,numCVReps); 
        vP2P2CVDist = NaN(1,numCVReps); vP2P1CVDist = NaN(1,numCVReps);
        vP1RCVDist = NaN(1,numCVReps); vP2RCVDist = NaN(1,numCVReps);
        
        % Loop
        curInd = 1; %Index of projection vectors (e.g., vPP)
        for CVRep = 1:numCVReps
            CVRep
            % Divide into folds 
            [trajStruct1,trajStruct2] = splitDataforCVMaxTrials(trajStruct,binWidth);  
            trajStruct1_P1 = trajStruct1([trajStruct1.posture]==1 | [trajStruct1.posture]==3);
            trajStruct1_P2 = trajStruct1([trajStruct1.posture]==4 | [trajStruct1.posture]==5);
            trajStruct2_P1 = trajStruct2([trajStruct2.posture]==1 | [trajStruct2.posture]==3);
            trajStruct2_P2 = trajStruct2([trajStruct2.posture]==4 | [trajStruct2.posture]==5);
            
            % Draw a random subspace
            v1 = normrnd(0,1,numPCsToKeep,1);
            DR = orth(v1);
            
            % Get P1 and P2 sig
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
            
            % Cross-projection
            vP1P1 = NaN(1,2); vP1P2 = NaN(1,2); 
            vP2P2 = NaN(1,2); vP2P1 = NaN(1,2);
            vP1R = NaN(1,2); vP2R = NaN(1,2);
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
                vP1R(projGroup) = trace(DR'*CP1*DR)/trace(CP1);
                vP2R(projGroup) = trace(DR'*CP2*DR)/trace(CP2);
            end
            
            vP1P1CVDist(CVRep) = mean(vP1P1);
            vP1P2CVDist(CVRep) = mean(vP1P2);
            vP2P2CVDist(CVRep) = mean(vP2P2);
            vP2P1CVDist(CVRep) = mean(vP2P1);
            vP1RCVDist(CVRep) = mean(vP1R);
            vP2RCVDist(CVRep) = mean(vP2R);
        end
        
        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset; 
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).vP1P1CVDist = vP1P1CVDist; 
        resultStruct(structInd).vP1P2CVDist = vP1P2CVDist; 
        resultStruct(structInd).vP2P2CVDist = vP2P2CVDist; 
        resultStruct(structInd).vP2P1CVDist = vP2P1CVDist; 
        resultStruct(structInd).vP1RCVDist = vP1RCVDist;
        resultStruct(structInd).vP2RCVDist = vP2RCVDist;
        
        resultStruct(structInd).vP1P1Mean = mean(vP1P1CVDist); 
        resultStruct(structInd).vP1P2Mean = mean(vP1P2CVDist); 
        resultStruct(structInd).vP2P1Mean = mean(vP2P1CVDist); 
        resultStruct(structInd).vP2P2Mean = mean(vP2P2CVDist); 
        
        structInd = structInd + 1; 
        
    end
    
%% Get monkey list 
    %Across datasets
    resultMonkeyList = {};
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
%% For each monkey, get overall distributions, CI's, and p-values
    monkeyResultStruct = struct('monkey',[],...
        'vP1P1CVMean',[],'vP1P1CI',[],...
        'vP1P2CVMean',[],'vP1P2CI',[],...
        'vP1RMean',[],'vP1RCI',[],...
        'vP1pVal',[],'vP1P2RpVal',[],...
        'vP2P2CVMean',[],'vP2P2CI',[],...
        'vP2P1CVMean',[],'vP2P1CI',[],...
        'vP2RMean',[],'vP2RCI',[],...
        'vP2pVal',[],'vP2P1RpVal',[]);
    
    monkeyResultStructInd = 1;
    for monkey = monkeyList
        %Get monkey data
        tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));        
        %Combine across datasets to get means and distributions
        vP1P1CVDist = [tempResultStruct.vP1P1CVDist]*100;
        vP1P1CVMean = mean(vP1P1CVDist);
        vP1P2CVDist = [tempResultStruct.vP1P2CVDist]*100;
        vP1P2CVMean = mean(vP1P2CVDist);
        vP1RCVDist = [tempResultStruct.vP1RCVDist]*100;
        
        vP2P2CVDist = [tempResultStruct.vP2P2CVDist]*100;
        vP2P2CVMean = mean(vP2P2CVDist);
        vP2P1CVDist = [tempResultStruct.vP2P1CVDist]*100;
        vP2P1CVMean = mean(vP2P1CVDist);
        vP2RCVDist = [tempResultStruct.vP2RCVDist]*100;
        
        %Test if vP1P1 significatly different from vP1P2 (paired, two-tailed)
        vP1PairedDiffDist = vP1P1CVDist-vP1P2CVDist;
        pd = fitdist(vP1PairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vP1pVal = 2*min([tail_1,tail_2]);
        
        %Test if vP1P2 significantly different from random (paired, two-tailed)
        vP1RPairedDiffDist = vP1P2CVDist-vP1RCVDist;
        pd = fitdist(vP1RPairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vP1P2RpVal = 2*min([tail_1,tail_2]);


        %Test if vP2P2 significantly different from vP2P1 (paired, two-tailed)
        vP2PairedDiffDist = vP2P2CVDist-vP2P1CVDist;
        pd = fitdist(vP2PairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vP2pVal = 2*min([tail_1,tail_2]);
        
        %Test if vP2P1 significantly different from random (unpaired, two-tailed)
        vP2RPairedDiffDist = vP2P1CVDist-vP2RCVDist;
        pd = fitdist(vP2RPairedDiffDist','Normal');
        tail_1 = cdf('Normal',0,pd.mu,pd.sigma);
        tail_2 = 1- cdf('Normal',0,pd.mu,pd.sigma);
        vP2P1RpVal = 2*min([tail_1,tail_2]);        
        
        
        %Add to monkeyResultStruct
        monkeyResultStruct(monkeyResultStructInd).monkey = monkey{1,1};
        
        monkeyResultStruct(monkeyResultStructInd).vP1P1CVMean = vP1P1CVMean;
        monkeyResultStruct(monkeyResultStructInd).vP1P1CI = [vP1P1CVMean-prctile(vP1P1CVDist,alpha*100/2) prctile(vP1P1CVDist,100-(alpha*100/2))-vP1P1CVMean];        
        monkeyResultStruct(monkeyResultStructInd).vP1P2CVMean = vP1P2CVMean;
        monkeyResultStruct(monkeyResultStructInd).vP1P2CI = [vP1P2CVMean-prctile(vP1P2CVDist,alpha*100/2) prctile(vP1P2CVDist,100-(alpha*100/2))-vP1P2CVMean];
        monkeyResultStruct(monkeyResultStructInd).vP1pVal = vP1pVal;
        monkeyResultStruct(monkeyResultStructInd).vP1RMean = mean(vP1RCVDist);
        monkeyResultStruct(monkeyResultStructInd).vP1RCI = [prctile(vP1RCVDist,alpha*100/2) prctile(vP1RCVDist,100-(alpha*100/2))]; 
        monkeyResultStruct(monkeyResultStructInd).vP1P2RpVal = vP1P2RpVal;
        
        monkeyResultStruct(monkeyResultStructInd).vP2P2CVMean = vP2P2CVMean;       
        monkeyResultStruct(monkeyResultStructInd).vP2P2CI = [vP2P2CVMean-prctile(vP2P2CVDist,alpha*100/2) prctile(vP2P2CVDist,100-(alpha*100/2))-vP2P2CVMean];        
        monkeyResultStruct(monkeyResultStructInd).vP2P1CVMean =vP2P1CVMean;      
        monkeyResultStruct(monkeyResultStructInd).vP2P1CI = [vP2P1CVMean-prctile(vP2P1CVDist,alpha*100/2) prctile(vP2P1CVDist,100-(alpha*100/2))-vP2P1CVMean];        
        monkeyResultStruct(monkeyResultStructInd).vP2pVal = vP2pVal;
        monkeyResultStruct(monkeyResultStructInd).vP2RMean = mean(vP2RCVDist);
        monkeyResultStruct(monkeyResultStructInd).vP2RCI = [prctile(vP2RCVDist,alpha*100/2) prctile(vP2RCVDist,100-(alpha*100/2))];       
        monkeyResultStruct(monkeyResultStructInd).vP2P1RpVal = vP2P1RpVal;
             

        
        monkeyResultStructInd = monkeyResultStructInd + 1;
    end
    
    %Export monkeyResultStruct
    if saveFig
        save(fullfile(saveDir,'monkeyResultStructS5.mat'),'monkeyResultStruct');
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
       vP1P1CVMean = tempPlotStruct.vP1P1CVMean;
       vP1P1CI = tempPlotStruct.vP1P1CI;
       vP1P2CVMean = tempPlotStruct.vP1P2CVMean;
       vP1P2CI = tempPlotStruct.vP1P2CI;
       
       bar(monkeyInd,vP1P1CVMean,'FaceColor',grey,'EdgeColor',black);
       hold on;
       errorbar(monkeyInd,vP1P1CVMean,vP1P1CI(1),vP1P1CI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0);      
       bar(monkeyInd+offset,vP1P2CVMean,'FaceColor',white,'EdgeColor',black);
       errorbar(monkeyInd+offset,vP1P2CVMean,vP1P2CI(1),vP1P2CI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'Capsize',0);
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
        saveas(gcf,fullfile(saveDir,'_P1SigVarExpl.svg'));
    end
    
    %Bar plot for P2 signal
    f = figure;  hold on;
    f.Position = [200 200 figWidth figHeight];
    monkeyInd = 1;
    for monkey = monkeyList
       tempPlotStruct = monkeyResultStruct(strcmp(plotMonkeyList,monkey{1,1}));
       vP2P2CVMean = tempPlotStruct.vP2P2CVMean;
       vP2P2CI = tempPlotStruct.vP2P2CI;
       vP2P1CVMean = tempPlotStruct.vP2P1CVMean;
       vP2P1CI = tempPlotStruct.vP2P1CI;
       
       bar(monkeyInd,vP2P1CVMean,'FaceColor',grey,'EdgeColor',black);
       errorbar(monkeyInd,vP2P1CVMean,vP2P1CI(1),vP2P1CI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0);
       hold on;
       
       bar(monkeyInd+offset,vP2P2CVMean,'FaceColor',white,'EdgeColor',black);
       errorbar(monkeyInd+offset,vP2P2CVMean,vP2P2CI(1),vP2P2CI(2),'Color',mcmap(monkeyInd,:),'LineWidth',3,'CapSize',0 );

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
        saveas(gcf,fullfile(saveDir,'P2SigVarExpl.svg'));
    end  