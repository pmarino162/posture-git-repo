clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\S3 - PC traces';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    %pcmap = vertcat(pcmap,pcmap(3:4,:));
    
%% Run Loop    
    tempList = {'N20171215','N20190226','R20201020','R20200221'};
    for datasetList = {'E20200316'}%{'E20200317','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}%{'R20200221'}%{'E20200317','E20200116','E20210706'}
       clf; close all
        %{'E20200317','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}
        % Load Data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318'}
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',-50},{'state','Step 1','first',250}};
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
            case 'E20210901'
                taskID = [Data.conditionData]; taskID = [taskID.taskID];
                Data = Data(taskID==1);
                trialInclStates(1).trialName = {'BCI Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
            %Iso
            case 'E20200116'
                trialInclStates(1).trialName = {'IsometricForce_1D'};   
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',500}};
            %Reaching
            case 'E20210706'
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
            case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
            case{'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
        end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        
        %Keep some postures for earl reaching
        switch dataset
            case {'E20210706'}
                    trajStruct = trajStruct(ismember([trajStruct.posture],[1,2,3,4,5,7]));
        end
        
        % Get minimum number of timestamps in condition averages
        numTimestamps = [];
        for i = 1:size(trajStruct,2)
            numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
        end
        [minNumTimestamps,i] = min(numTimestamps);
               
        switch dataset
            case {'E20200316'}
                pcmap = flip(pcmap);
            case {'N20171215'}     
                pcmap = vertcat(pcmap(3,:),pcmap(1,:),pcmap(5,:));
            case {'E20210901'}
                pcmap = vertcat(pcmap(1,:),pcmap(3,:),pcmap(3,:),pcmap(5,:),pcmap(4,:));
                pcmap = vertcat(pcmap(4,:),pcmap(5,:),pcmap(5,:),pcmap(1,:),pcmap(3,:));
            case {'R20201020'}
                pcmap = vertcat(pcmap(2,:),pcmap(1,:));
        end
        
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
        
        %Do PCA on all condition averages
        allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
        [allPCA,score,latent,tsquared,explained,allMu] = pca(allTraj);
        totalVar = trace(cov(allTraj));
        
        %Add Projections to trajStruct
        for i = 1:size(trajStruct,2)
            %Add average traces to trajStruct
            trajStruct(i).PCA.traj = (trajStruct(i).avgSmoothFR.traj-allMu)*allPCA;
            trajStruct(i).PCA.VAF =  100.*(diag(cov(allTraj*allPCA))')./totalVar;
            trajStruct(i).PCA.time = trajStruct(i).avgSmoothFR.timestamps;
            %Compute CI, add to trajStruct
            %adjust CI so that it uses available number of trials at each time point
            numTrials = size(trajStruct(i).allSmoothFR,2);
            maxNumTimestamps = 0;
            for trial = 1:numTrials
               numTimestamps = size(trajStruct(i).allSmoothFR(trial).traj,1);
               if numTimestamps > maxNumTimestamps
                   maxNumTimestamps = numTimestamps;
               end
            end
            tempAllTraj = NaN(maxNumTimestamps,numTrials,numChannels);
            for trial = 1:numTrials
                traj = trajStruct(i).allSmoothFR(trial).traj;
                numTimestamps = size(traj,1);
                tempAllTraj(1:numTimestamps,trial,:) = traj*allPCA;
            end
            for timestamp = 1:maxNumTimestamps
                timestampMat = squeeze(tempAllTraj(timestamp,:,:)); %All trials for all channels for current timestamp
                tempNumTrials = sum(~isnan(timestampMat(:,1)));
                trajStruct(i).PCA.CI(timestamp,:) = 1.96.*nanstd(timestampMat)./sqrt(tempNumTrials);
            end
        end
        
        for plotGroup = 1:2
            %Set up plotting params
            if plotGroup == 1                           
                plotPostureList = [1,5];
                plotTargetList = [1,5];
            else
                plotPostureList = [1,5];
                plotTargetList = [1];
            end
            if strcmpi(dataset,'E20210706')
                plotPostureList = [2,7];
            end
            if strcmpi(dataset,'N20171215') || strcmpi(dataset,'N20190226') || strcmpi(dataset,'R20201020') || strcmpi(dataset,'R20200221')
                plotPostureList = [1,2];
            end

            %Plot traces
            fs = 14;
            f = figure; f.Position = [50 50 1000 400];
            for posture = plotPostureList
                for target = plotTargetList
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.traj;
                    CI = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.CI;
                    time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.time;                    
                    for dim = 1:6
                        subplot(2,3,dim); hold on;
                        if target == 1
                            shadedErrorBar(time,traj(:,dim),CI(1:length(time),dim),'lineprops',{'LineWidth',2','Color',pcmap(posture,:)});
                        elseif target == 5
                            shadedErrorBar(time,traj(:,dim),CI(1:length(time),dim),'lineprops',{'--','LineWidth',2','Color',pcmap(posture,:)});
                        end
                    end
                end
            end
            maxYRange = 0;
            for dim = 1:6
                subplot(2,3,dim)
                ax = gca;
                ylimits = ax.YLim;
                yRange = ylimits(2)-ylimits(1);
                if yRange > maxYRange
                    maxYRange = yRange;
                end
            end
            for dim = 1:6
                VAF = round(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.VAF);
                subplot(2,3,dim);
                ax = gca;
                ylimits = ax.YLim;
                yMid = ylimits(1) + 0.5*(ylimits(2)-ylimits(1));
                ylim([yMid-maxYRange/2 yMid+maxYRange/2])
                ylabel(['PC ',num2str(dim),' (',num2str(VAF(dim)),'%)'])
                yticks([min(ax.YTick) max(ax.YTick)]);
                if dim < 4
                    xticks([])
                else
                    xticks([min(ax.XTick) max(ax.XTick)]);
                end
                if dim == 5
                    xlabel('time (ms)')
                end
                set(gca,'fontname','arial')
                set(gca,'fontsize',fs)
            end           
            if saveFig
                if plotGroup == 1
                    saveas(gcf,fullfile(saveDir,dataset,'PC_traces_2targets.svg'));
                else
                    saveas(gcf,fullfile(saveDir,dataset,'PC_traces_1target.svg'));
                end
            end

            %Plot - Orthographic
            fs = 14;
            figure
            %f = figure; f.Position = [50 50 200 200];
            hold on
            timePts = 1:minNumTimestamps;
            xDim = 2; yDim = 3; zDim = 1;
            for posture = plotPostureList
                for target = plotTargetList
                    if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.traj(timePts,:); 
                        if target == 1
                            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                        elseif target == 5
                            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'--','Color',pcmap(posture,:),'LineWidth',2);
                        end

                       % plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    end
                end
            end


            grid on
            xticklabels({}); yticklabels({}); zticklabels({}); 
            VAF = round(trajStruct(1).PCA.VAF);

            xlabel(['PC',num2str(xDim),' (',num2str(VAF(xDim)),'%)'])
            ylabel(['PC',num2str(yDim),' (',num2str(VAF(yDim)),'%)'])
            zlabel(['PC',num2str(zDim),' (',num2str(VAF(zDim)),'%)'])
            view([20 10])
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
            if saveFig
                if plotGroup == 1
                    saveas(gcf,fullfile(saveDir,dataset,'PC_Orth_2targets.fig'))
                    saveas(gcf,fullfile(saveDir,dataset,'PC_Orth_2targets.svg'));
                else
                    saveas(gcf,fullfile(saveDir,dataset,'PC_Orth_1target.fig'));
                    saveas(gcf,fullfile(saveDir,dataset,'PC_Orth_1target.svg'));
                end
            end
        end
    end