clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
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
                structInd = structInd + 1;
            end
        end
        
        %Plot and save
            fs = 14;
            %Posture marginalizations for blocks 1 and 2 in block 1
            %posture PC space
            postureList = unique([trajStruct.posture]);
            tempPcmap = pcmap;
            if length(postureList) == 3
                tempPcmap = pcmap([1,3,5],:);
            end
            %Individual trial observations in block 1 posture PC space
            figure; hold on;
            for i = 1:numel(trialStruct)
               trial = trialStruct(i).trial;
               posture = trialStruct(i).posture;
               PC = trialStruct(i).PC;
               plot(trial,PC,'.','MarkerSize',10,'Color',tempPcmap(posture,:));
            end
            %Add posture means for E20190830 dataset
            if strcmpi(dataset,'E20190830')
                for posture = 1:5
                   allPts = [trialStruct([trialStruct.posture]==posture).PC]; 
                   avg = mean(allPts);
                   plot(1400,avg,'.','MarkerSize',30,'Color',tempPcmap(posture,:));
                end
                xticks([0 600 1200])
            end
            xlabel('trial'); ylabel('Posture Dim 1 (a.u.)')
            %title(dataset)
            set(gca, 'TickDir', 'out')
            set(gca,'fontname','arial'); set(gca,'fontsize',fs)
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_postureAxisVTrial.svg']));
            end
            

            
    end
    
%% Local function 
function [postureMargTraj,posturePCs] = getMarginalsAndPCs(trajStruct)
   
    %Get minimum number of timestamps in condition averages
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    [minNumTimestamps,i] = min(numTimestamps);

    %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);

    %Form X, containing trial-averaged data for each condition
    X = NaN(minNumTimestamps,numTargets,numPostures,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            X(:,targetInd,postureInd,:) = traj(1:minNumTimestamps,:); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end

    % Do marginalizations of X
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;

    %Compute PCs of marginalizations
    targetMargTraj = squeeze(targetMargTraj);
        targetMargTraj = reshape(targetMargTraj,[numTargets*minNumTimestamps,numChannels]);
        [targetPCs,~,~,~,~,targetMargMu] = pca(targetMargTraj); 
    postureMargTraj = squeeze(postureMargTraj);
        postureMargTraj = reshape(postureMargTraj,[numPostures*minNumTimestamps,numChannels]);
        [posturePCs,~,~,~,~,postureMargMu] = pca(postureMargTraj); 
        
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
    postureMargTraj = squeeze(postureMargTraj);
end

