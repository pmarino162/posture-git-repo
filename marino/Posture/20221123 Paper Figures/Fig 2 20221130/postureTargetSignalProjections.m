clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'D:\Posture Paper\20221123\Figures\Figure 2\Subspace Projections';
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
    resultStructInd = 1;
    for datasetList = {'E20200116'}%{'E20200116'}%'E20210901'}%,'E20200317','E20200318'}%,'N20171215','N20180221','R20201020','R20201021'}%,'R20201020','R20201021'}    %{'E20200317','N20180221'}      
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
        xticklabels({}); yticklabels({});
        grid on; set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_TSigTSub.svg']));
        end
                
        %Plot P signal in P subspace
        figure
        pDim1 = 1; pDim2 = 2;
        if strcmp(dataset,'E20210901')
            pDim1 = 2; pDim2 = 3;
        end
        hold on
        postureInd = 1;
        for posture = postureList
            traj = postureMargTraj((postureInd-1)*minNumTimestamps+1:postureInd*minNumTimestamps,:);
            postureProj = traj*posturePCs;
            plot(-postureProj(:,pDim1),postureProj(:,pDim2),'LineWidth',2,'Color',pcmap(posture,:));
            postureInd = postureInd + 1;
        end
        xlabel(['Posture Dim ',num2str(pDim1)]); ylabel(['Posture Dim ',num2str(pDim2)])
        grid on; xticklabels({}); yticklabels({});
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PSigPSub.svg']));
        end

        %Plot T and P signal in T subspace 
        figure; hold on
        targetInd = 1;
        for target = targetList
            traj = targetMargTraj((targetInd-1)*minNumTimestamps+1:targetInd*minNumTimestamps,:);
            targetProj = traj*targetPCs;
            plot(targetProj(:,1),targetProj(:,2),'LineWidth',2,'Color','k');
            targetInd = targetInd +1;
        end
        postureInd = 1;
        for posture = postureList
            traj = postureMargTraj((postureInd-1)*minNumTimestamps+1:postureInd*minNumTimestamps,:);
            postureProj = traj*targetPCs;
            plot(postureProj(:,pDim1),postureProj(:,pDim2),'LineWidth',2,'Color',pcmap(posture,:));
            postureInd = postureInd + 1;
        end
        xlabel('Target Dim 1'); ylabel('Target Dim 2')
        xticklabels({}); yticklabels({});
        grid on; set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_TPSigTSub.svg']));
        end
        
        %Plot T and P signal in P subspace 
        figure; hold on
        targetInd = 1;
        for target = targetList
            traj = targetMargTraj((targetInd-1)*minNumTimestamps+1:targetInd*minNumTimestamps,:);
            targetProj = traj*posturePCs;
            plot(targetProj(:,1),targetProj(:,2),'LineWidth',2,'Color','k');
            targetInd = targetInd + 1;
        end
        postureInd = 1;  
        for posture = postureList
            traj = postureMargTraj((postureInd-1)*minNumTimestamps+1:postureInd*minNumTimestamps,:);
            postureProj = traj*posturePCs;
            plot(-postureProj(:,pDim1),postureProj(:,pDim2),'LineWidth',2,'Color',pcmap(posture,:));
            postureInd = postureInd + 1;
        end
        xlabel(['Posture Dim ',num2str(pDim1)]); ylabel(['Posture Dim ',num2str(pDim2)])
        grid on; xticklabels({}); yticklabels({});
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_TPSigPSub.svg']));
        end
        
        %Plot full trajectories in joint subspace
        figure; hold on;
        postureInd = 1; 
        for posture = postureList
            for target = [1,5]
                targetInd = find(targetList==target);
                traj = squeeze(X(:,targetInd,postureInd,:));
                traj = traj-CIMargTraj;
                traj = traj*PTOrth;
                plot3(traj(:,1),traj(:,2),traj(:,3),'LineWidth',2,'Color',pcmap(posture,:));
                plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
            postureInd = postureInd + 1;
        end
        axis equal; grid on;
        xlabel('Posture Dim 1'); ylabel('Posture Dim 2'); zlabel('Target Dim 1')
        set(gca,'FontSize',fs); set(gca,'FontName','Arial');
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PTOrth.svg']));
            saveas(gcf,fullfile(saveDir,[dataset,'_PTOrth.fig']));
        end
        

        % Plot cross-projection eigenspectra   
        fs = 20;
        figure
            b = bar([targetTPC_VAF(1:numTargetPCs),postureTPC_VAF(1:numTargetPCs)],'FaceColor','flat');
            for k = 1:2
               b(k).CData = .7*(k-1)*[1 1 1]; 
            end
            xlabel('Target PC'); ylabel('Variance Explained (%)')
            legend('Target Signal','Posture Signal')
            set(gca, 'TickDir', 'out'); set(gca,'fontname','arial'); set(gca,'fontsize',fs)
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_TargetPCSpec.svg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_TargetPCSpec.fig']));
            end
        figure
            b = bar([targetPPC_VAF(1:numPosturePCs),posturePPC_VAF(1:numPosturePCs)],'FaceColor','flat');
            for k = 1:2
               b(k).CData = .7*(k-1)*[1 1 1]; 
            end
            xlabel('Posture PC'); ylabel('Variance Explained (%)')
            set(gca, 'TickDir', 'out'); set(gca,'fontname','arial'); set(gca,'fontsize',fs)
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_PosturePCSpec.svg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_PosturePCSpec.fig']));
            end

       
    end
