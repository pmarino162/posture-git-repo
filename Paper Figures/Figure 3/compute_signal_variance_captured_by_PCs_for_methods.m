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
        'explainedP_2dims',[],'explainedT_2dims',[]);
    structInd = 1;
    numTotalPostures = 0;
    numRemainingPostures = 0;
    %Run loop
    for datasetList = [reachDatasetList,bciDatasetList,isoDatasetList]
        %% Set up trajStruct
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);            
        [postureList,~,targetList,~,~,~] = getTrajStructDimensions(trajStruct);
        numTotalPostures = numTotalPostures + length(postureList);
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
        numRemainingPostures = numRemainingPostures + length(postureList);
 
        %% Split data into two groups for cross-validation
        [minNumCondTrials] = getMinNumCondTrials(trajStruct);
        [trajStruct1,trajStruct2] = splitDataforCV(trajStruct,minNumCondTrials,binWidth);  
        
        %% Fit P and T spaces for each group, compute projection variances
        [minNumTimestamps] = getMinNumTimestamps(trajStruct);    
        [pSig1,tSig1] = getPandTsig(trajStruct1,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
            pSig1Reshape = reshape(squeeze(pSig1),[numPostures*minNumTimestamps,numPCsToKeep]);
            tSig1Reshape = reshape(squeeze(tSig1),[numTargets*minNumTimestamps,numPCsToKeep]);
            [pDimsGroup1,~,~,~,explainedP_group1,pSigMu] = pca(pSig1Reshape); 
            [tDimsGroup1,~,~,~,explainedT_group1,tSigMu] = pca(tSig1Reshape); 
        [pSig2,tSig2] = getPandTsig(trajStruct2,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
            pSig2Reshape = reshape(squeeze(pSig2),[numPostures*minNumTimestamps,numPCsToKeep]);
            tSig2Reshape = reshape(squeeze(tSig2),[numTargets*minNumTimestamps,numPCsToKeep]);
            [pDimsGroup2,~,~,~,explainedP_group2,~] = pca(pSig2Reshape); 
            [tDimsGroup2,~,~,~,explainedT_group2,~] = pca(tSig2Reshape); 
        
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).explainedP_2dims = [sum(explainedP_group1(1:2)), sum(explainedP_group2(1:2))];
        resultStruct(structInd).explainedT_2dims = [sum(explainedT_group1(1:2)), sum(explainedT_group2(1:2))];
        structInd = structInd + 1;
        
    end
      
%% Compute overall mean and SD    
allExplainedP = [];
allExplainedT = [];
for i = 1:size(resultStruct,2)
    allExplainedP = [allExplainedP,resultStruct(i).explainedP_2dims];
    allExplainedT = [allExplainedT,resultStruct(i).explainedT_2dims];
end
mean(allExplainedP)
std(allExplainedP)
mean(allExplainedT)
std(allExplainedT)

%% Compute fraction of postures remaining after exclusion
(numRemainingPostures/numTotalPostures)*100