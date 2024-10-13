clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\Reviews\Analyses\condition_trial_counts';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Datasets to include in analysis 
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};         
    isoDatasetList = {'E20200116','E20200117','E20200120'};    

%% Set up result struct
    resultStruct = struct('monkey',[],'dataset',[],'numCondTrials',[]);
    structInd = 1;

%% Main loop
    task = 'Reach';
    if strcmp(task,'BCI')
        taskDatasetList = bciDatasetList;
    elseif strcmp(task,'Iso')
        taskDatasetList = isoDatasetList;
    elseif strcmp(task,'Reach')
        taskDatasetList = reachDatasetList;
    end

    for datasetList = taskDatasetList%{'E20200316'}%{'N20190222','N20190226','R20200221','R20200222'}%reachDatasetList%{'E20200316'}%bciDatasetList% reachDatasetList%{'E20210707','N20190226','R20200221'}%{'E20200316','N20171215','R20201020'}%{ 'R20200221'}%bciDatasetList%{'E20210706','E20210707','E20210708','E20210709'}%reachDatasetList
        %% Set up trajStruct
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',false);  
        %Remove any conditions for which there weren't enough trials
        [numCondTrials] = getNumCondTrials(trajStruct,'showHist',false);    
        %Fill result struct
        resultStruct(structInd).monkey = dataset(1);
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).numCondTrials = numCondTrials;
        structInd = structInd + 1;
    end
%% Get monkey list
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).monkey;
    end
    monkeyList = unique(resultMonkeyList);
    
    
%% Make histogram
for monkey = monkeyList
    figure
    hold on;
    
    tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
    
    for i = 1:size(tempResultStruct,2)
        histogram([tempResultStruct(i).numCondTrials],'BinWidth',2,'DisplayStyle','stairs')
    end
    legend(tempResultStruct.dataset)
    title(monkey{1,1})
    xlabel('Num. Trials')
    ylabel('Num. Conditions')
    saveas(gcf,fullfile(saveDir,[monkey{1,1},'_',task,'_condition_trial_counts.svg'])); 
end