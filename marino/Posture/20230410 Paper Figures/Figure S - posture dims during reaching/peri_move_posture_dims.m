clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20230203\orthographic dpca pop traj';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Run Loop    
    datasetList = {'N20190226'}%{'E20210706'}%{'N20190226',
       clf; close all
        % Load Data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            case{'N20190226'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
            case{'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};

        end
                        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
                periMoveTrajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-100},{'kin','moveOnsetTime','first',300}};
                allTrajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
                trialInclStates(1).inclStates = {{'state','Target Hold','first',0},{'state','Success with Reward','first',0}};
                targetHoldTrajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        trajStruct = periMoveTrajStruct;
        
        %Keep postures 1-4 only for earl
        switch dataset
            case {'E20210706'}
                    trajStruct = trajStruct(ismember([trajStruct.posture],[1,2,3,4]));
        end
        

        
        % Get minimum number of timestamps in condition averages
        numTimestamps = [];
        for i = 1:size(trajStruct,2)
            numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
        end
        [minNumTimestamps,i] = min(numTimestamps);
        
        

        
        
%         %Custom timestamps/pcmap
%         if minNumTimestamps > 6
%             minNumTimestamps = 6;
%         end
        
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
        

        Xdpca = permute(X,[4,3,2,1]);

        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        margNames = {'P', 'T', 'CI', 'PTI'};
        margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

        N = numChannels;
        P = numPostures;
        T = numTargets;
        timePts = minNumTimestamps;
        time = trajStruct(1).avgSmoothFR.timestamps(1:minNumTimestamps);
        timeEvents = [];

        [W,V,whichMarg] = dpca(Xdpca, 10, ...
            'combinedParams', combinedParams);

        explVar = dpca_explainedVariance(Xdpca, W, V, ...
            'combinedParams', combinedParams);

        dpca_plot(Xdpca, W, V, @dpca_plot_default, ...
            'explainedVar', explVar, ...
            'marginalizationNames', margNames, ...
            'marginalizationColours', margColours, ...
            'whichMarg', whichMarg,                 ...
            'time', time,                        ...
            'timeEvents', timeEvents,               ...
            'timeMarginalization', 3, ...
            'legendSubplot', 1);


        % Plot Marginal dims
        postureMargID = find(strcmp(margNames,'P'));
        targetMargID = find(strcmp(margNames,'T'));
        CIMargID = find(strcmp(margNames,'CI'));

        postureDims = find(whichMarg==postureMargID,2);
        targetDims = find(whichMarg==targetMargID,2);
        CIDims = find(whichMarg==CIMargID,2);
        %dimList = [postureDims(1),targetDims(1),CIDims(1),postureDims(2),targetDims(2),CIDims(2)];

        plotTargetList = [1,5];
        if strcmpi(dataset,'E20200116')
            plotTargetList = [3,7];
        end
        
        
        %Orthonormalize pdpca1 and tdpca1-2
        pDPCA = W(:,postureDims);
        tDPCA = W(:,targetDims);
        [PTOrth,~] = qr([pDPCA(:,1),tDPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
        
        
        
        %Add Projections to trajStruct
        totalVar = trace(cov(allTraj));
        for i = 1:size(periMoveTrajStruct,2)
            %Add average traces to trajStruct
            periMoveTrajStruct(i).PTOrth.traj = periMoveTrajStruct(i).avgSmoothFR.traj*PTOrth;
            periMoveTrajStruct(i).tDPCA.traj = periMoveTrajStruct(i).avgSmoothFR.traj*tDPCA;
            periMoveTrajStruct(i).pDPCA.traj = periMoveTrajStruct(i).avgSmoothFR.traj*pDPCA;
            %Get VAF
            periMoveTrajStruct(i).PTOrth.VAF =  [explVar.componentVar(postureDims(1)),explVar.componentVar(targetDims(1)),explVar.componentVar(targetDims(2))];%100.*(diag(cov(allTraj*PTOrth))')./totalVar;
        end

        for i = 1:size(allTrajStruct,2)
            %Add average traces to trajStruct
            allTrajStruct(i).PTOrth.traj = allTrajStruct(i).avgSmoothFR.traj*PTOrth;
            allTrajStruct(i).tDPCA.traj = allTrajStruct(i).avgSmoothFR.traj*tDPCA;
            allTrajStruct(i).pDPCA.traj = allTrajStruct(i).avgSmoothFR.traj*pDPCA;
            %Get VAF
            allTrajStruct(i).PTOrth.VAF =  [explVar.componentVar(postureDims(1)),explVar.componentVar(targetDims(1)),explVar.componentVar(targetDims(2))];%100.*(diag(cov(allTraj*PTOrth))')./totalVar;
        end
        
       for i = 1:size(targetHoldTrajStruct,2)
            %Add average traces to trajStruct
            targetHoldTrajStruct(i).PTOrth.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*PTOrth;
            targetHoldTrajStruct(i).tDPCA.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*tDPCA;
            targetHoldTrajStruct(i).pDPCA.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*pDPCA;
            %Get VAF
            targetHoldTrajStruct(i).PTOrth.VAF =  [explVar.componentVar(postureDims(1)),explVar.componentVar(targetDims(1)),explVar.componentVar(targetDims(2))];%100.*(diag(cov(allTraj*PTOrth))')./totalVar;
        end
        
        
        %Plot - Orthographic
        trajStruct = periMoveTrajStruct;
        fs = 14;
        figure
        hold on
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1;
        for posture = postureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                    plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                end
            end
        end
             
        
        grid on
        xticklabels({}); yticklabels({}); zticklabels({}); 
        VAF = round(trajStruct(1).PTOrth.VAF);

        xlabel(['Goal Dim 1 (',num2str(VAF(xDim)),'%)'])
        ylabel(['Goal Dim 2 (',num2str(VAF(yDim)),'%)']) 
        zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
        view([20 10])
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PTOrth.fig']));
        end
        
        
        plotTargetList = [4,8];
        trajStruct = allTrajStruct;
        figure; hold on;
        for posture = postureList
            for target = plotTargetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).pDPCA.traj(:,1); 
                    time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps; 
                    plot(time,traj,'Color',pcmap(posture,:),'LineWidth',2);
                end
            end
        end
        
                plotTargetList = [4,8];
        trajStruct = targetHoldTrajStruct;
        figure; hold on;
        for posture = postureList
            for target = plotTargetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).pDPCA.traj(:,1); 
                    time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps; 
                    plot(time,traj,'Color',pcmap(posture,:),'LineWidth',2);
                end
            end
        end