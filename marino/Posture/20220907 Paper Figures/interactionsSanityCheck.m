clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'D:\Plots\Interactions Sanity Check';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = customRainbow;
     
%% Run loop for each dataset   
    resultStructInd = 1;
    for datasetList = {'N20190222','N20190226','N20190227','N20190228','R20200221','R20200222'}%,'N20171215','N20180221','R20201020','R20201021'}%,'R20201020','R20201021'}    %{'E20200317','N20180221'}      
        close all
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            case {'E20200316','E20200317','E20200318'}
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Release','first',0},{'state','Target Hold','first',0}};
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
            %Earl Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
            %Nigel and Rocky Reaching
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305',...
                    'N20190306','N20190307','R20200221','R20200222'}
                if strcmpi(dataset(1),'N')
                    trialInclStates(1).trialName = {'Nigel Dissociation'};
                elseif strcmpi(dataset(1),'R')
                    trialInclStates(1).trialName = {'Rocky Dissociation'};
                end
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
            end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
            
        % Get minimum number of timestamps in condition averages
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
        
        %Do PCA on all trajectories
        allTraj = [];
        for i = 1:size(trajStruct,2)
            traj = trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps,:);
            allTraj = vertcat(allTraj,traj);
        end
        [allPCs,~,~,~,explained,allMu] = pca(allTraj);
        
        %Project into PCs
        for i = 1:size(trajStruct,2)
            traj = trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps,:);
            trajStruct(i).PCA.traj = (traj-allMu)*allPCs(:,1:10);
        end
        
        switch dataset
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                postureList = 1:4;
        end
        %Plot and save figures
        for target = targetList
            maxYLim = -100; minYLim = 100;
            f = figure;
            f.Position = [100 100 1000 500];
            for posture = postureList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.traj; 
                    time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(1:minNumTimestamps); 
                    for dim = 1:10
                        subplot(2,5,dim)
                        plot(time,traj(:,dim),'Color',pcmap(posture,:),'LineWidth',2);
                        hold on
                        ax = gca;
                        ylim = ax.YLim;
                        if ylim(1) < minYLim 
                            minYLim = ylim(1);
                        end
                        if ylim(2) > maxYLim 
                            maxYLim = ylim(2);
                        end
                        if posture == postureList(end)
                           ax.YLim = [minYLim maxYLim];
                           xlabel('time (ms)')
                           ylabel(['PC ',num2str(dim)])
                        end
                    end
                end
            end
            sgtitle(['Target ',num2str(target)])
            saveas(gcf,fullfile(saveDir,dataset,['T',num2str(target),'_PCs.svg']));
        end
        
        
    end
        