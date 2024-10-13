clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 3';
    %saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S5 - different joints';
    set(0, 'DefaultFigureRenderer', 'painters');
     
%% Set Parameters
   numPCsToKeep = 10;     %Num PCs to project data into before analysis
   
%% Datasets to include in analysis 
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','E20210901','N20171215','N20180221','R20201020','R20201021'};         
    isoDatasetList = {'E20200116','E20200117','E20200120'};
    
%% Main loop   
    resultStruct = struct('animal',[],'dataset',[],'VAF',[]);
    structInd = 1;
    for curDataset = [reachDatasetList,bciDatasetList,isoDatasetList]
        %% Set up trajStruct
        %Load data
        dataset = curDataset{1,1}
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);            
        %Project all data down to top PCs
        [allTraj] = collectAllAvgTraj(trajStruct);
        [trajStruct] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep);
        
        VAF = sum(trajStruct(1).avgPCA.VAF)
        
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).VAF = VAF;
        structInd = structInd + 1;
    end
%% Compute overall mean and sd
allVAF = [];
for i = 1:size(resultStruct,2)
    allVAF = [allVAF,resultStruct(i).VAF];
end
mean(allVAF)
std(allVAF)
