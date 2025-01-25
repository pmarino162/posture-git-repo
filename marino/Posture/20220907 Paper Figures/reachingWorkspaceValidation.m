clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Plots\Interactions Sanity Check';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    tcmap = customRainbow;
     
%% Run loop for each dataset   
    resultStructInd = 1;
    for datasetList = {'N20190226'}   
        
        close all
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR','marker','markerVel'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            %Earl Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',200}};
            %Nigel and Rocky Reaching
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305',...
                    'N20190306','N20190307','R20200221','R20200222'}
                if strcmpi(dataset(1),'N')
                    trialInclStates(1).trialName = {'Nigel Dissociation'};
                elseif strcmpi(dataset(1),'R')
                    trialInclStates(1).trialName = {'Rocky Dissociation'};
                end
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',200}};
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
        
        %Plot - all marker traj in workspace
        figure; hold on;
        for posture = postureList
            for target = targetList
                numCondTrials = size(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allMarker,2);
                for trial = 1:numCondTrials
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allMarker(trial).traj; 
                    plot(traj(:,1),traj(:,2),'Color',pcmap(posture,:),'LineWidth',1);
                end
            end
        end
        xlabel('x (mm)'); ylabel('y (mm)');
        axis equal
        
        %Plot - marker vel in workspace (avg)
        figure; hold on;
        for posture = postureList
            for target = targetList
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarkerVel.traj; 
                plot(traj(:,1),traj(:,2),'Color',pcmap(posture,:),'LineWidth',1);
            end
        end
        xlabel('Vx'); ylabel('Vy');
        axis equal
        
        %Plot - pos and vel vs. time for each posture (all)
        for posture = postureList
            figure; 
            for target = targetList
                numCondTrials = size(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allMarker,2);
                for trial = 1:numCondTrials
                    pos = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allMarker(trial).traj; 
                    vel = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allMarkerVel(trial).traj; 
                    time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allMarker(trial).timestamps; 
                    subplot(4,1,1)
                        plot(time,pos(:,1),'Color',tcmap(target,:),'LineWidth',1); hold on;
                        ylabel('x')
                    subplot(4,1,2)
                        plot(time,pos(:,2),'Color',tcmap(target,:),'LineWidth',1); hold on;
                        ylabel('y')
                    subplot(4,1,3)
                        plot(time,vel(:,1),'Color',tcmap(target,:),'LineWidth',1); hold on;
                        ylabel('vx')
                    subplot(4,1,4)
                        plot(time,vel(:,2),'Color',tcmap(target,:),'LineWidth',1); hold on;
                        ylabel('vy')
                        xlabel('time rel. MO (ms)')
                end
            end
        end
        
        
        %Plot - pos and vel vs. time for each posture (avg)
        for posture = postureList
            figure; 
            for target = targetList
                pos = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj; 
                vel = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarkerVel.traj; 
                time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps; 
                subplot(4,1,1)
                    plot(time,pos(:,1),'Color',tcmap(target,:),'LineWidth',1); hold on;
                    ylabel('x')
                subplot(4,1,2)
                    plot(time,pos(:,2),'Color',tcmap(target,:),'LineWidth',1); hold on;
                    ylabel('y')
                subplot(4,1,3)
                    plot(time,vel(:,1),'Color',tcmap(target,:),'LineWidth',1); hold on;
                    ylabel('vx')
                subplot(4,1,4)
                    plot(time,vel(:,2),'Color',tcmap(target,:),'LineWidth',1); hold on;
                    ylabel('vy')
                    xlabel('time rel. MO (ms)')
            end
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
            %saveas(gcf,fullfile(saveDir,dataset,['T',num2str(target),'_PCs.svg']));
        end
        
        
    end
      