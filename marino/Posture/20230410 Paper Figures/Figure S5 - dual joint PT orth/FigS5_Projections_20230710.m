clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\S5 - dual-joint PT orth';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],'Psig_to_Tspace_AI',[],'Psig_to_Tspace_nullDist',[],'Tsig_to_Pspace_AI',[],'Tsig_to_Pspace_nullDist',[]);
    
%% Parameters
    % How many PC's to describe the 'manifold'? 
    manVarThreshold = 90;
    % Keep Enough PCs to capture 90% of respective variance
    varThreshold = 90;
    
%% Run loop for each dataset      
    %Load data
    dataset = 'E20210901';
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
        case {'E20200116'} 
            trialInclStates(1).trialName = {'IsometricForce_1D'};   
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',200}};
        case {'E20210706'}
            trialInclStates(1).trialName = {'GridReaching'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
        case 'E20210901'
            taskID = [Data.conditionData]; taskID = [taskID.taskID];
            Data = Data(taskID==1);
            trialInclStates(1).trialName = {'BCI Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
        case {'N20171215','N20180221'}
            trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
        case {'R20201020','R20201021'}
            trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

    %Remove postures for earl reach
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
        
    %Do PCA; mean center; project down to manifold
    allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
    [allPCs,~,~,~,explained,allMu] = pca(allTraj);
    numManPCs = 1;
    while sum(explained(1:numManPCs)) < manVarThreshold
        numManPCs = numManPCs + 1;
    end
    X = NaN(minNumTimestamps,numTargets,numPostures,numManPCs);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            X(:,targetInd,postureInd,:) = (traj(1:minNumTimestamps,:)-allMu)*allPCs(:,1:numManPCs); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end

    % Do marginalizations of X
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Offsets
    targetMargOffset = mean(X,[1 3],'omitnan') - CIMargOffset;
    postureMargOffset = mean(X,[1 2],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
    %targetMargTraj = mean(X,[3],'omitnan');
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
    targetTrajNoCP = X - CIMargTraj - CIMargOffset - postureMargTraj - postureMargOffset;

    %Compute PCs of marginalizations
    CIMargTraj = squeeze(CIMargTraj);
        CIMargOffset = squeeze(CIMargOffset);
        [CIPCs,~,~,~,~,CIMu] = pca(CIMargTraj); 
    targetMargTraj = squeeze(targetMargTraj);
        targetMargTraj = reshape(targetMargTraj,[numTargets*minNumTimestamps,numManPCs]);
        [targetPCs,~,~,~,~,targetMargMu] = pca(targetMargTraj); 
    targetTrajNoCP = squeeze(targetTrajNoCP);
        targetTrajNoCP = reshape(targetTrajNoCP,[numTargets*numPostures*minNumTimestamps,numManPCs]);
        [targetNoCPPCs,~,~,~,~,targetNoCPMargMu] = pca(targetTrajNoCP); 
    postureMargTraj = squeeze(postureMargTraj);
        postureMargTraj = reshape(postureMargTraj,[numPostures*minNumTimestamps,numManPCs]);
        [posturePCs,~,~,~,~,postureMargMu] = pca(postureMargTraj); 

    % Get total and cross-projection VAF for each marginalization
    totalTargetVar = trace(cov(targetMargTraj));
    totalPostureVar = trace(cov(postureMargTraj));
    targetTPC_VAF = diag(cov(targetMargTraj*targetPCs))./totalTargetVar.*100;
    postureTPC_VAF = diag(cov(postureMargTraj*targetPCs))./totalPostureVar.*100;
    targetPPC_VAF = diag(cov(targetMargTraj*posturePCs))./totalTargetVar.*100;
    posturePPC_VAF = diag(cov(postureMargTraj*posturePCs))./totalPostureVar.*100;
        
    %Get number of each type of PC to include for subsequent analysis
    numTargetPCs = 1;
        while sum(targetTPC_VAF(1:numTargetPCs)) < varThreshold
            numTargetPCs = numTargetPCs + 1;
        end
    numPosturePCs = 1;
        while sum(posturePPC_VAF(1:numPosturePCs)) < varThreshold
            numPosturePCs = numPosturePCs + 1;
        end
    
       %Form joint subspace for full traj visualization      
       [PTOrth,~] = qr([posturePCs(:,1:2),targetPCs(:,1)]); PTOrth = PTOrth(:,1:3);    
            
              
        fs = 14;
        
        %Plot T signal in T subspace    
        figure; hold on
        targetInd = 1;
        for target = targetList
            traj = targetMargTraj((targetInd-1)*minNumTimestamps+1:targetInd*minNumTimestamps,:);
            targetProj = traj*targetPCs;
            plot(targetProj(:,1),targetProj(:,2),'LineWidth',2,'Color','k');
            targetInd = targetInd + 1;
        end
        xlabel('Target Dim 1'); ylabel('Target Dim 2')
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        axis equal; 
        %grid on; xticklabels({}); yticklabels({});
        set(gca,'TickDir','out');
        ax = gca; TSubXLims = ax.XLim; TSubYLims = ax.YLim;
        TSubXRange = TSubXLims(2)-TSubXLims(1); TSubYRange = TSubYLims(2)-TSubYLims(1);
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_TSigTSub.svg']));
        end
                     
        
        %Plot P signal in P subspace
        figure; hold on
        postureInd = 1;
        for posture = postureList
            traj = postureMargTraj((postureInd-1)*minNumTimestamps+1:postureInd*minNumTimestamps,:);
            postureProj = traj*posturePCs;
            plot(-postureProj(:,1),postureProj(:,2),'LineWidth',2,'Color',pcmap(posture,:));
            postureInd = postureInd + 1;
        end
        xlabel('Posture Dim 1'); ylabel('Posture Dim 2')
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        axis equal; 
                %grid on; xticklabels({}); yticklabels({});
        set(gca,'TickDir','out');
        ax = gca; PSubXLims = ax.XLim; PSubYLims = ax.YLim;
        PSubXRange = PSubXLims(2)-PSubXLims(1); PSubYRange = PSubYLims(2)-PSubYLims(1);
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PSigPSub.svg']));
        end

        
        %Plot 3D Psignal
        figure; hold on
        postureInd = 1;
        for posture = postureList
            traj = postureMargTraj((postureInd-1)*minNumTimestamps+1:postureInd*minNumTimestamps,:);
            postureProj = traj*posturePCs;
            plot3(-postureProj(:,1),postureProj(:,2),postureProj(:,3),'LineWidth',2,'Color',pcmap(posture,:));
            postureInd = postureInd + 1;
        end
        xlabel('Posture Dim 1'); 
        ylabel('Posture Dim 2'); 
        zlabel('Posture Dim 3');
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        axis equal; 
                %grid on; xticklabels({}); yticklabels({});
        set(gca,'TickDir','out');
        view([-25,21])
        grid on
        xticklabels({}); yticklabels({}); zticklabels({}); 
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PSigPSub3D.svg']));
        end
        
        
        
        
        
        %Plot P signal in T subspace 
        figure; hold on
        postureInd = 1;
        for posture = postureList
            traj = postureMargTraj((postureInd-1)*minNumTimestamps+1:postureInd*minNumTimestamps,:);
            postureProj = traj*targetPCs;
            plot(postureProj(:,1),postureProj(:,2),'LineWidth',2,'Color',pcmap(posture,:));
            postureInd = postureInd + 1;
        end
        xlabel('Target Dim 1'); ylabel('Target Dim 2')
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        ax = gca; xlim = ax.XLim; ylim = ax.YLim;
        xmid = mean(xlim); ymid = mean(ylim);
        ax.XLim = [xmid-PSubXRange/2 xmid+PSubXRange/2]; 
        ax.YLim = [ymid-PSubYRange/2 ymid+PSubYRange/2]; 
                %grid on; xticklabels({}); yticklabels({});
        set(gca,'TickDir','out');
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PSigTSub.svg']));
        end

        
        %Plot T signal in P subspace 
        figure; hold on
        targetInd = 1;
        for target = targetList
            traj = targetMargTraj((targetInd-1)*minNumTimestamps+1:targetInd*minNumTimestamps,:);
            targetProj = traj*posturePCs;
            plot(targetProj(:,1),targetProj(:,2),'LineWidth',2,'Color','k');
            targetInd = targetInd + 1;
        end
        xlabel('Posture Dim 1'); ylabel('Posture Dim 2')
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        ax = gca; xlim = ax.XLim; ylim = ax.YLim;
        xmid = mean(xlim); ymid = mean(ylim);
        ax.XLim = [xmid-TSubXRange/2 xmid+TSubXRange/2]; 
        ax.YLim = [ymid-TSubYRange/2 ymid+TSubYRange/2]; 
                %grid on; xticklabels({}); yticklabels({});
        set(gca,'TickDir','out');
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_TSigPSub.svg']));
        end
        

       
