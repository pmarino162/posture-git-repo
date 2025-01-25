clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\S3 - PC traces';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Get trajStruct
    dataset = 'E20200316';
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
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

      
    % Get minimum number of timestamps in condition averages
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    [minNumTimestamps,i] = min(numTimestamps);

    switch dataset
        case {'E20200316'}
            pcmap = flip(pcmap);
    end

    %Get posture and target lists
    postureList = [1,3,5];% unique([trajStruct.posture]);
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
        

%% Plot traces
    plotPostureList = postureList;
    plotTargetList = [1:8];%1:8;
    fs = 4;
    
    for target = plotTargetList
        f = figure; f.Position = [50 50 200 80];
        for posture = plotPostureList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.traj;
            CI = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.CI;
            time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.time;  
            for dim = 1:6
                subplot(2,3,dim); hold on;
                shadedErrorBar(time,traj(:,dim),CI(1:length(time),dim),'lineprops',{'LineWidth',1,'Color',pcmap(posture,:)});
            end
        end    
        
        %Format
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
            ax.TickDir = 'out';
            if dim < 4
                xticks([])
            else
                xticks([-50 0 250]);
            end
            if dim == 5
                %xlabel('time (ms)')
            end
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
        end           

        
        %Save
        if saveFig
            saveas(gcf,fullfile(saveDir,dataset,['target_',num2str(target),'_PC_traces.svg']));
        end
    end
    

    


            
