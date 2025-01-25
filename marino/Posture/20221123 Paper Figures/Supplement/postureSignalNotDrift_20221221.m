clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'D:\Posture Paper\20221123\Figures\Supplement\Posture Signal Is Not Slow Drift';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Run loop for each dataset    
    for datasetList = {'N20171215','E20190830'}%,'E20190830'}%{'N20180221','N20171215'}
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            case 'E20190830'
                trialInclStates(1).trialName = {'CenterOut_ForceBar_BC'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','BC Freeze','first',0},{'state','Success with Reward','first',0}};
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
        end

        %Get trajStruct for each block
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

        %Get posturePCs on blocks 1 and 2
        [postureMargTraj,posturePCs] = getMarginalsAndPCs(trajStruct);

        %Get postureLDA on block1
        [posturePCLDA] = getPostureBlockPCLDA(trajStruct);
        
        %Create trialStruct with individual trial observations
        trialStruct = struct('trial',[],'posture',[],'PC',[],'PCLDA',[]);
        structInd = 1;
        for i = 1:size(trajStruct,2)
            posture = trajStruct(i).posture;
            for j = 1:size(trajStruct(i).allSmoothFR,2)
                trial = trajStruct(i).allSmoothFR(j).trialNum;
                traj = trajStruct(i).allSmoothFR(j).traj;
                traj = mean(traj,1);
                trialStruct(structInd).trial = trial;
                trialStruct(structInd).posture = posture;
                trialStruct(structInd).PC = traj*posturePCs(:,1);
                trialStruct(structInd).PCLDA = traj*posturePCLDA(:,1);
                structInd = structInd + 1;
            end
        end
        
        %Plot and save
            fs = 14;
            %Posture marginalizations for blocks 1 and 2 in block 1
            %posture PC space
            postureList = unique([trajStruct.posture]);
            
            %Individual trial observations in block 1 posture PC space
            figure; hold on;
            for i = 1:numel(trialStruct)
               trial = trialStruct(i).trial;
               posture = trialStruct(i).posture;
               PC = trialStruct(i).PC;
               plot(trial,PC,'.','MarkerSize',10,'Color',pcmap(posture,:));
            end
            xlabel('trial'); ylabel('Posture PC 1')
            title(dataset)
            set(gca,'fontname','arial'); set(gca,'fontsize',fs)
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_postureAxisVTrial.svg']));
            end
            
            
            %Individual trial observations in block 1 posture PCLDA space
            figure; hold on;
            for i = 1:numel(trialStruct)
               trial = trialStruct(i).trial;
               posture = trialStruct(i).posture;
               PCLDA = trialStruct(i).PCLDA;
               plot(trial,PCLDA,'.','MarkerSize',10,'Color',pcmap(posture,:));
            end
            xlabel('trial'); ylabel('Posture PCLDA 1')
            
    end
    

    
