clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    %saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\BCIPopTraj';
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Figure 2';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Run loop for each dataset    
    for datasetList = {'E20200317'}%{'N20180221','N20171215'}
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            case 'E20200317'
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
        end
        
        %Get block labels
        blockLabels = [Data.conditionData];
        blockLabels = [blockLabels.block];
        numBlocks = max(blockLabels);
        
        %Get trajStruct for each block
        allTrajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        
        %Get posturePCs on blocks 1 and 2
        [block1PostureMargTraj,block1PosturePCs] = getMarginalsAndPCs(allTrajStruct);

        %Get postureLDA on block1
        [block1PosturePCLDA] = getPostureBlockPCLDA(allTrajStruct);
        
        %Create trialStruct with individual trial observations
        trialStruct = struct('trial',[],'posture',[],'b1PC',[],'b1PCLDA',[]);
        structInd = 1;
        for i = 1:size(allTrajStruct,2)
            posture = allTrajStruct(i).posture;
            for j = 1:size(allTrajStruct(i).allSmoothFR,2)
                trial = allTrajStruct(i).allSmoothFR(j).trialNum;
                traj = allTrajStruct(i).allSmoothFR(j).traj;
                traj = mean(traj,1);
                trialStruct(structInd).trial = trial;
                trialStruct(structInd).posture = posture;
                trialStruct(structInd).b1PC = traj*block1PosturePCs(:,1);
                trialStruct(structInd).b1PCLDA = traj*block1PosturePCLDA(:,1);
                structInd = structInd + 1;
            end
        end
        
        
        
        %Plot and save
            fs = 14;
            
            %Individual trial observations in block 1 posture PC space
            figure; hold on;
            for i = 1:numel(trialStruct)
               trial = trialStruct(i).trial;
               posture = trialStruct(i).posture;
               b1PC = trialStruct(i).b1PC;
               plot(trial,b1PC,'.','MarkerSize',10,'Color',pcmap(posture,:));
            end
            xlabel('trial')
            ylabel('Posture PC 1')
            
            %Individual trial observations in block 1 posture PCLDA space
            figure; hold on;
            for i = 1:numel(trialStruct)
               trial = trialStruct(i).trial;
               posture = trialStruct(i).posture;
               b1PCLDA = trialStruct(i).b1PCLDA;
               plot(trial,b1PCLDA,'.','MarkerSize',10,'Color',pcmap(posture,:));
            end
            xlabel('trial')
            ylabel('Posture PCLDA 1')
            
    end
    

    
