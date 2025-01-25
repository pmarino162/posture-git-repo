clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\ReachOrthQuant';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;

%% Run loop for each dataset    
    holdEpoch = 'center';
    for datasetList = {'E20210706','N20190226'} 
        %datasetList = {'N20190226'}
        %datasetList = {'E20210706'} 
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

%         %Look at distribution of center hold times; make holdData trials
%         %bigger than 25th pctile
%         figure
%         kinData = [Data.kinData];
%         centerHoldTimes = [kinData.centerHoldTime];
%         histogram(centerHoldTimes);
%         centerHoldTimeThreshold = prctile(centerHoldTimes,25);
%         holdData = Data(centerHoldTimes > centerHoldTimeThreshold);
%         xlabel('Center Hold Time (ms)')
%         ylabel('Count')
        
%         %Look at distribution of target hold times; set threshold
%         figure
%         kinData = [Data.kinData];
%         targetHoldTimes = [kinData.targetHoldTime];
%         histogram(targetHoldTimes);
%         xlabel('Target Hold Time (ms)')
%         ylabel('Count')
%         targetHoldTimeThreshold = prctile(targetHoldTimes,25);
%         %targetHoldTimeThreshold = 800;
%         targetHoldData = Data(targetHoldTimes > targetHoldTimeThreshold);

        centerHoldData = Data;
        targetHoldData = Data;
        
        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            case 'E20210706'
                %All
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',500}};
                %Reach
                reachTrialInclStates(1).trialName = {'GridReaching'};
                reachCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                reachTrialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',-50}};
                %Center Hold
                centerHoldTrialInclStates(1).trialName = {'GridReaching'};
                centerHoldCondFields = {{'posture','conditionData','postureID'}};
                centerHoldTrialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Delay','first',50}};
                %Target Hold
                targetHoldTrialInclStates(1).trialName = {'GridReaching'};
                targetHoldCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                targetHoldTrialInclStates(1).inclStates = {{'state','Target Hold','first',0},{'state','Success with Reward','first',0}};
            case 'N20190226'
                %All
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',500}};
                %Reach
                reachTrialInclStates(1).trialName = {'Nigel Dissociation'};
                reachCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                reachTrialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',-50}};
                %Center Hold
                centerHoldTrialInclStates(1).trialName = {'Nigel Dissociation'};
                centerHoldCondFields = {{'posture','conditionData','postureID'}};
                centerHoldTrialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Reach','first',50}};
                %Target Hold
                targetHoldTrialInclStates(1).trialName = {'Nigel Dissociation'};
                targetHoldCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                targetHoldTrialInclStates(1).inclStates = {{'state','Target Hold','first',0},{'state','Success with Reward','first',0}};

        end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        reachTrajStruct = getTrajStruct20220419(Data,reachCondFields,trajFields,reachTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        centerHoldTrajStruct = getTrajStruct20220419(centerHoldData,centerHoldCondFields,trajFields,centerHoldTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        targetHoldTrajStruct = getTrajStruct20220419(targetHoldData,targetHoldCondFields,trajFields,targetHoldTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        
        %Clear unneccessary variables
        clearvars Data centerHoldData targetHoldData
        % Get minimum number of timestamps in condition averages
        %All 
            numTimestamps = [];
            for i = 1:size(trajStruct,2)
                numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
            end
            figure
            histogram(numTimestamps)
            xlabel('Number of 25ms bins')
            ylabel('Number of conditions')
            [minNumTimestamps,i] = min(numTimestamps);
        %Reach
            reachNumTimestamps = [];
            for i = 1:size(reachTrajStruct,2)
                reachNumTimestamps(i) = length(reachTrajStruct(i).avgSmoothFR.timestamps);
            end
            figure
            histogram(reachNumTimestamps)
            xlabel('Number of 25ms bins')
            ylabel('Number of conditions')
            [reachMinNumTimestamps,i] = min(reachNumTimestamps);
        %Center Hold 
            centerHoldNumTimestamps = [];
            for i = 1:size(centerHoldTrajStruct,2)
                centerHoldNumTimestamps(i) = length(centerHoldTrajStruct(i).avgSmoothFR.timestamps);
            end
            figure
            histogram(centerHoldNumTimestamps)
            xlabel('Number of 25ms bins')
            ylabel('Number of conditions')
            [centerHoldMinNumTimestamps,i] = min(centerHoldNumTimestamps);
        %Target Hold
            targetHoldNumTimestamps = [];
            for i = 1:size(targetHoldTrajStruct,2)
                targetHoldNumTimestamps(i) = length(targetHoldTrajStruct(i).avgSmoothFR.timestamps);
            end
            figure
            histogram(targetHoldNumTimestamps)
            xlabel('Number of 25ms bins')
            ylabel('Number of conditions')
            [targetHoldMinNumTimestamps,i] = min(targetHoldNumTimestamps);
            
         %Get posture and target lists 
           %All 
            postureList = unique([trajStruct.posture]);
            numPostures = size(postureList,2);
            targetList = unique([trajStruct.target]); 
            numTargets = size(targetList,2);
            %Reach
            reachPostureList = unique([reachTrajStruct.posture]);
            reachNumPostures = size(reachPostureList,2);
            reachTargetList = unique([reachTrajStruct.target]); 
            reachNumTargets = size(reachTargetList,2);
            %Center Hold
            centerHoldPostureList = unique([centerHoldTrajStruct.posture]);
            centerHoldNumPostures = size(centerHoldPostureList,2); 
            %Target Hold
            targetHoldPostureList = unique([targetHoldTrajStruct.posture]);
            targetHoldNumPostures = size(targetHoldPostureList,2);
            targetHoldTargetList = unique([targetHoldTrajStruct.target]); 
            targetHoldNumTargets = size(targetHoldTargetList,2);
            
        %num Channels and numConditions    
            numChannels = size(reachTrajStruct(1).avgSmoothFR.traj,2);
            numConditions = size(reachTrajStruct,2);
        
        %Remove postures for which not all targets were attempted
            rmPostures = [];
            for posture = reachPostureList
                tempTrajStruct = reachTrajStruct([reachTrajStruct.posture]==posture);
                postureTargetList = [tempTrajStruct.target];
                if ~isequal(postureTargetList,reachTargetList)
                    rmPostures = [rmPostures,posture];
                end
            end
            if ~isempty(rmPostures)
                %All
                rmInd = [trajStruct.posture]==rmPostures;
                trajStruct(rmInd) = [];
                postureList(postureList==rmPostures) = [];
                numPostures = length(postureList);
                %Reach
                rmInd = [reachTrajStruct.posture]==rmPostures;
                reachTrajStruct(rmInd) = [];
                reachPostureList(reachPostureList==rmPostures) = [];
                reachNumPostures = length(reachPostureList);
                %Target Hold
                rmInd = [targetHoldTrajStruct.posture]==rmPostures;
                targetHoldTrajStruct(rmInd) = [];
                targetHoldPostureList(targetHoldPostureList==rmPostures) = [];
                targetHoldNumPostures = length(targetHoldPostureList);              
            end
          
        
        %X
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
        
        %reachX
        reachX = NaN(reachMinNumTimestamps,reachNumTargets,reachNumPostures,numChannels);
        postureInd = 1;
        for posture = reachPostureList
            targetInd = 1;
            for target = reachTargetList
                traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgSmoothFR.traj;
                reachX(:,targetInd,postureInd,:) = traj(1:reachMinNumTimestamps,:); 
                targetInd = targetInd + 1;
            end
            postureInd = postureInd+1;
        end
        
        %centerHoldX
        centerHoldX = NaN(centerHoldMinNumTimestamps,centerHoldNumPostures,numChannels);
        postureInd = 1;
        for i=1:size(centerHoldTrajStruct,2)
            traj = centerHoldTrajStruct(i).avgSmoothFR.traj;
            centerHoldX(:,postureInd,:) = traj(1:centerHoldMinNumTimestamps,:); 
            postureInd = postureInd+1;
        end
        
        %targetHoldX
        targetHoldX = NaN(targetHoldMinNumTimestamps,targetHoldNumTargets,targetHoldNumPostures,numChannels);
        postureInd = 1;
        for posture = targetHoldPostureList
            targetInd = 1;
            for target = targetHoldTargetList
                traj = targetHoldTrajStruct([targetHoldTrajStruct.posture]==posture & [targetHoldTrajStruct.target]==target).avgSmoothFR.traj;
                targetHoldX(:,targetInd,postureInd,:) = traj(end-targetHoldMinNumTimestamps+1:end,:); 
                targetInd = targetInd + 1;
            end
            postureInd = postureInd+1;
        end
             
        % Do marginalizations of reachX
        %reachXdims: 1=time, 2=target, 3=posture, 4=channel
        %Condition-invariant
        reachCIMargOffset = mean(reachX,[1 2 3],'omitnan');
        reachCIMargTraj = mean(reachX,[2 3],'omitnan') - reachCIMargOffset;
        %Posture and Target Offsets
        reachTargetMargOffset = mean(reachX,[1 3],'omitnan') - reachCIMargOffset;
        reachPostureMargOffset = mean(reachX,[1 2],'omitnan') - reachCIMargOffset;
        %Posture and Target Traj
        reachTargetMargTraj = mean(reachX,[3],'omitnan') - reachCIMargTraj  - reachCIMargOffset;
        reachPostureMargTraj = mean(reachX,[2],'omitnan') - reachCIMargTraj - reachCIMargOffset;
        reachTargetTrajNoCP = reachX - reachCIMargTraj - reachCIMargOffset - reachPostureMargTraj - reachPostureMargOffset;
        reachNoCI = reachX - reachCIMargTraj - reachCIMargOffset;
        
        % Do marginalizations of targetHoldX
        %targetHoldXdims: 1=time, 2=target, 3=posture, 4=channel
        %Condition-invariant
        targetHoldCIMargOffset = mean(targetHoldX,[1 2 3],'omitnan');
        targetHoldCIMargTraj = mean(targetHoldX,[2 3],'omitnan') - targetHoldCIMargOffset;
        %%Posture and Target Offsets
        %targetHoldTargetMargOffset = mean(targetHoldX,[1 3],'omitnan') - targetHoldCIMargOffset;
        %targetHoldPostureMargOffset = mean(targetHoldX,[1 2],'omitnan') - targetHoldCIMargOffset;
        %%Posture and Target Traj
        %targetHoldTargetMargTraj = mean(targetHoldX,[3],'omitnan') - targetHoldCIMargTraj  - targetHoldCIMargOffset;
        %targetHoldPostureMargTraj = mean(targetHoldX,[2],'omitnan') - targetHoldCIMargTraj - targetHoldCIMargOffset;
        %targetHoldTargetTrajNoCP = targetHoldX - targetHoldCIMargTraj - targetHoldCIMargOffset - targetHoldPostureMargTraj - targetHoldPostureMargOffset;
        targetHoldNoCI = targetHoldX - targetHoldCIMargTraj - targetHoldCIMargOffset;
        
        % Do marginalizations of centerHoldX
        %centerHoldXdims: 1=time, 2=posture, 3=channel
        %Condition-invariant
        centerHoldCIMargOffset = mean(centerHoldX,[1 2],'omitnan');
        centerHoldCIMargTraj = mean(centerHoldX,[2],'omitnan') - centerHoldCIMargOffset;
        %Posture Offset
        centerHoldPostureMargOffset = mean(centerHoldX,[1 2],'omitnan') - centerHoldCIMargOffset;
        %Posture Traj
        centerHoldNoCI = centerHoldX - centerHoldCIMargTraj - centerHoldCIMargOffset; 
        
       
        
        %Do PCA on all condition averaged trajectories
        allReachTraj = reshape(reachX,[reachNumTargets*reachNumPostures*reachMinNumTimestamps,numChannels]);
        allTargetHoldTraj = reshape(targetHoldX,[targetHoldNumTargets*targetHoldNumPostures*targetHoldMinNumTimestamps,numChannels]);
        allCenterHoldTraj = reshape(centerHoldX,[centerHoldNumPostures*centerHoldMinNumTimestamps,numChannels]);
        
        if strcmp(holdEpoch,'center')
            allTraj = vertcat(allReachTraj,allCenterHoldTraj);
        elseif strcmp(holdEpoch,'target')
            allTraj = vertcat(allReachTraj,allTargetHoldTraj);
        end
        
        [allPCA,score,latent,tsquared,explained,mu] = pca(allTraj); 
        [reachPCA,score,latent,tsquared,explained,mu] = pca(allReachTraj); 
        [centerHoldPCA,score,latent,tsquared,explained,mu] = pca(allCenterHoldTraj); 
        [targetHoldPCA,score,latent,tsquared,explained,mu] = pca(allTargetHoldTraj); 
        Ca = cov(allTraj);
        [Ua,Sa] = eig(Ca);
 

        %Define CIMargTraj, postureMargTraj, targetMargTraj, targetTrajNoCP
        CIMargTraj = reachCIMargTraj;
        if strcmp(holdEpoch,'center')
            postureMargTraj = centerHoldNoCI;
        elseif strcmp(holdEpoch,'target')
            postureMargTraj = targetHoldNoCI;
        end
        targetMargTraj = reachTargetMargTraj;
        targetTrajNoCP = reachTargetTrajNoCP;
        
        CIMargTraj = squeeze(CIMargTraj);
        [CIPCA,score,latent,tsquared,explained,mu] = pca(CIMargTraj); 

        targetMargTraj = squeeze(targetMargTraj);
        targetMargTraj = reshape(targetMargTraj,[reachNumTargets*reachMinNumTimestamps,numChannels]);
        [Dt,score,latent,tsquared,explained,mu] = pca(targetMargTraj); 

        targetTrajNoCP = squeeze(targetTrajNoCP);
        targetTrajNoCP = reshape(targetTrajNoCP,[reachNumTargets*reachNumPostures*reachMinNumTimestamps,numChannels]);
        [targetNoCPPCA,score,latent,tsquared,explained,mu] = pca(targetTrajNoCP); 

        postureMargTraj = squeeze(postureMargTraj);
        if strcmp(holdEpoch,'center')
            postureMargTraj = reshape(postureMargTraj,[centerHoldNumPostures*centerHoldMinNumTimestamps,numChannels]);
        elseif strcmp(holdEpoch,'target')
            postureMargTraj = reshape(postureMargTraj,[targetHoldNumTargets*targetHoldNumPostures*targetHoldMinNumTimestamps,numChannels]);
        end
        [Dp,score,latent,tsquared,explained,mu] = pca(postureMargTraj); 

        % Get cross-projection VAF
        totalTargetVar = trace(cov(targetMargTraj));
        totalPostureVar = trace(cov(postureMargTraj));
        targetTPC_VAF = diag(cov(targetMargTraj*Dt))./totalTargetVar.*100;
        postureTPC_VAF = diag(cov(postureMargTraj*Dt))./totalPostureVar.*100;
        targetPPC_VAF = diag(cov(targetMargTraj*Dp))./totalTargetVar.*100;
        posturePPC_VAF = diag(cov(postureMargTraj*Dp))./totalPostureVar.*100;
        
        % Keep Enough PCs to capture 90% of respective variance
        varThreshold = 90;

        numTargetPCs = 1;
        while sum(targetTPC_VAF(1:numTargetPCs)) < varThreshold
            numTargetPCs = numTargetPCs + 1;
        end

        numPosturePCs = 1;
        while sum(posturePPC_VAF(1:numPosturePCs)) < varThreshold
            numPosturePCs = numPosturePCs + 1;
        end
    
        % Visualize full trajectories in PT PC subspace
        figure
        if strcmp(holdEpoch,'center')
            pcmap = parula(centerHoldNumPostures);
        elseif strcmp(holdEpoch,'target')
            pcmap = parula(targetHoldNumPostures*targetHoldNumTargets);
        end
        hold on
        dims = [1,2,2];
        postureInd = 1;
        for posture = reachPostureList
            for target = reachTargetList
                traj = squeeze(reachNoCI(:,target,postureInd,:));
                targetProj = traj*Dt;
                postureProj = traj*Dp;
                plot3(targetProj(:,dims(1)),targetProj(:,dims(2)),postureProj(:,dims(3)),'LineWidth',2,'Color',pcmap(posture,:));
            end
            postureInd = postureInd + 1;
        end
        
        
        %Visualize target marginalizations in target space 
        figure
        pcmap = parula(7);
        hold on
        dims = [1,2,3];
        for target = reachTargetList
                traj = squeeze(reachTargetMargTraj(:,target,:));
                targetTraj = traj*Dt;
                plot3(targetTraj(:,dims(1)),targetTraj(:,dims(2)),targetTraj(:,dims(3)),'LineWidth',2,'Color',pcmap(posture,:));
        end
        
        %Visualize posture marginalizations in Posture space
        figure
        hold on
        dims = [1,2,3];
        if strcmp(holdEpoch,'center')
            pcmap = parula(centerHoldNumPostures);
            for posture = centerHoldPostureList
                    traj = squeeze(centerHoldNoCI(:,posture,:));
                    postureProj = traj*Dp;
                    plot3(postureProj(:,dims(1)),postureProj(:,dims(2)),postureProj(:,dims(3)),'LineWidth',2,'Color',pcmap(posture,:));
            end
        elseif strcmp(holdEpoch,'target')
            postureInd = 1;
            for posture = targetHoldPostureList
                pcmap = parula(targetHoldNumPostures*targetHoldNumTargets);
                for target = targetHoldTargetList
                    traj = squeeze(targetHoldNoCI(:,target,postureInd,:));
                    postureProj = traj*Dp;
                    plot3(postureProj(:,dims(1)),postureProj(:,dims(2)),postureProj(:,dims(3)),'LineWidth',2,'Color',pcmap(posture,:));
                end
                postureInd = postureInd + 1;
            end
        end
        
        
%         %Visualize posture marginalizations in Target space
%         figure
%         pcmap = parula(7);
%         hold on
%         dims = [1,2,3];
%         for posture = centerHoldPostureList
%                 traj = squeeze(centerHoldNoCI(:,posture,:));
%                 targetProj = traj*Dt;
%                 plot3(targetProj(:,dims(1)),targetProj(:,dims(2)),targetProj(:,dims(3)),'LineWidth',2,'Color',pcmap(posture,:));
%         end
        
        
        
        
        
        % Plot cross-projection eigenspectra   
        figure
            bar([targetTPC_VAF(1:numTargetPCs),postureTPC_VAF(1:numTargetPCs)])
            ax = gca;
            ax.FontSize = 16;
            xlabel('Target PC')
            ylabel('Variance Explained (%)')
            legend('Target Marginalization','Posture Marginalization')
            set(gca, 'TickDir', 'out')
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_TargetPCSpec.svg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_TargetPCSpec.fig']));
            end
            
        figure
            if numPosturePCs == 1
                bar([targetPPC_VAF(1:2),posturePPC_VAF(1:2)])
                xlim([0.5 1.5])
            else
                bar([targetPPC_VAF(1:numPosturePCs),posturePPC_VAF(1:numPosturePCs)])
            end
            ax = gca;
            ax.FontSize = 16;
            xlabel('Posture PC')
            ylabel('Variance Explained (%)')
            set(gca, 'TickDir', 'out')  
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_PosturePCSpec.svg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_PosturePCSpec.fig']));
            end
        
            
        %Cumulative version of cross-projections
        fs = 14;
        cumTargetTPC_VAF = targetTPC_VAF;
        cumPostureTPC_VAF = postureTPC_VAF;
        for i = 2:numTargetPCs
           cumTargetTPC_VAF(i) = cumTargetTPC_VAF(i) + cumTargetTPC_VAF(i-1);
           cumPostureTPC_VAF(i) = cumPostureTPC_VAF(i) + cumPostureTPC_VAF(i-1);
        end
        cumTargetPPC_VAF = targetPPC_VAF;
        cumPosturePPC_VAF = posturePPC_VAF;
        for i = 2:numPosturePCs
           cumTargetPPC_VAF(i) = cumTargetPPC_VAF(i) + cumTargetPPC_VAF(i-1);
           cumPosturePPC_VAF(i) = cumPosturePPC_VAF(i) + cumPosturePPC_VAF(i-1);
        end
        figure
            bar([cumTargetTPC_VAF(1:numTargetPCs),cumPostureTPC_VAF(1:numTargetPCs)])
            hold on
                plot([0 numTargetPCs-.1],[cumTargetTPC_VAF(numTargetPCs) cumTargetTPC_VAF(numTargetPCs)],'--','LineWidth',2,'Color',[0, 0.4470, 0.7410])
                plot([0 numTargetPCs+.1],[cumPostureTPC_VAF(numTargetPCs) cumPostureTPC_VAF(numTargetPCs)],'--','LineWidth',2,'Color',[0.8500, 0.3250, 0.0980])
            yticks([0 cumPostureTPC_VAF(numTargetPCs) cumTargetTPC_VAF(numTargetPCs) 100])
            yticklabels({'0',num2str(round(cumPostureTPC_VAF(numTargetPCs))),num2str(round(cumTargetTPC_VAF(numTargetPCs))),'100'})
            xlabel('Number of Target Dimensions Included')
            ylabel('Cumulative Target Variance Explained (%)')
            legend('Target Signal','Posture Signal')
            set(gca, 'TickDir', 'out')
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_CumulativeTargetPCSpec.svg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_CumulativeTargetPCSpec.fig']));
            end
        figure
            bar([cumTargetPPC_VAF(1:numPosturePCs),cumPosturePPC_VAF(1:numPosturePCs)])
            hold on
                plot([0 numPosturePCs-.1],[cumTargetPPC_VAF(numPosturePCs) cumTargetPPC_VAF(numPosturePCs)],'--','LineWidth',2,'Color',[0, 0.4470, 0.7410])
                plot([0 numPosturePCs+.1],[cumPosturePPC_VAF(numPosturePCs) cumPosturePPC_VAF(numPosturePCs)],'--','LineWidth',2,'Color',[0.8500, 0.3250, 0.0980])
            yticks([0 cumTargetPPC_VAF(numPosturePCs) cumPosturePPC_VAF(numPosturePCs) 100])
            yticklabels({'0',num2str(round(cumTargetPPC_VAF(numPosturePCs))),num2str(round(cumPosturePPC_VAF(numPosturePCs))),'100'})
            xlabel('Number of Posture Dimensions Included')
            ylabel('Cumulative Posture Variance Explained (%)')
            set(gca, 'TickDir', 'out')  
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_CumulativePosturePCSpec.svg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_CumulativePosturePCSpec.fig']));
            end
            
            
        %How much target-related variance is captured in posture subspace?  What's the most you could hope to capture w an equal number of dims?
        targetPCaptured = sum(diag(cov(targetMargTraj*Dp(:,1:numPosturePCs))));
        targetTCaptured = sum(diag(cov(targetMargTraj*Dt(:,1:numPosturePCs))));
    
        %How much posture-related variance is captured in target subspace?  What's the most you could hope to capture w an equal number of dims?
%         if numTargetPCs > size(Dp,2)
%             numTargetPCs = size(Dp,2);
%         end
        postureTCaptured = sum(diag(cov(postureMargTraj*Dt(:,1:numTargetPCs))));
        posturePCaptured = sum(diag(cov(postureMargTraj*Dp(:,1:numTargetPCs))));

        %Alignment Index
            alignT = targetPCaptured/targetTCaptured;
            alignP = postureTCaptured/posturePCaptured;
    
        %Generate random subspaces; ask what fraction of target-related and posture-related variability is captured by them
            numDraws = 10000;

            %Target
            randTVarCaptured = zeros(1,numDraws);
            for i = 1:numDraws
               %Draw a random orthonormalized subspace that matches data covariance
               %& has #posturePCs dims
               v = normrnd(0,1,numChannels,numPosturePCs);
               Dr = orth((Ua*sqrt(Sa)*v)/norm(Ua*sqrt(Sa)*v));
               %Compute amount of target-related and posture-related variability that is captured by that
               %space
               randTVarCaptured(i) = sum(diag(cov(targetMargTraj*Dr)));
            end
            randTAlignmentDist = randTVarCaptured/targetTCaptured;

            %Posture
            randPVarCaptured = zeros(1,numDraws);
            for i = 1:numDraws
               %Draw a random orthonormalized subspace that matches data covariance
               %& has #targetPCs dims
               v = normrnd(0,1,numChannels,numTargetPCs);
               Dr = orth((Ua*sqrt(Sa)*v)/norm(Ua*sqrt(Sa)*v));
               %Compute amount of target-related and posture-related variability that is captured by that
               %space
               randPVarCaptured(i) = sum(diag(cov(postureMargTraj*Dr)));
            end
            randPAlignmentDist = randPVarCaptured/posturePCaptured;
    
        %Plot p-value histograms     
            figure
                hist(randTAlignmentDist)
                hold on
                ax = gca;
                xLim = ax.XLim; yLim = ax.YLim;
                line([alignT,alignT],yLim,'Color','r','LineWidth',2);
                pT = sum(randTAlignmentDist <= alignT)/numDraws;
                text(mean(xLim),yLim(2)-(yLim(2)-yLim(1))*.1,['p = ',num2str(pT)])
                xlabel('Alignment Index')
                ylabel('Count')

             figure
                hist(randPAlignmentDist)
                hold on
                ax = gca;
                xLim = ax.XLim; yLim = ax.YLim;
                line([alignP,alignP],yLim,'Color','r','LineWidth',2);
                pP = sum(randPAlignmentDist <= alignP)/numDraws;
                text(mean(xLim),yLim(2)-(yLim(2)-yLim(1))*.1,['p = ',num2str(pP)])
                xlabel('Alignment Index')
                ylabel('Count')

              randTAlignmentMean = mean(randTAlignmentDist);
              randTAlignmentCI = 1.96*std(randTAlignmentDist)/sqrt(numDraws);

              randPAlignmentMean = mean(randPAlignmentDist);
              randPAlignmentCI = 1.96*std(randPAlignmentDist)/sqrt(numDraws);

            %Plot results
                f = figure;
                f.Position = [100 100 500 500];
                bar([alignT,randTAlignmentMean,alignP,randPAlignmentMean])
                hold on
                errorbar(2,randTAlignmentMean,randTAlignmentCI,'k','LineWidth',2)
                errorbar(4,randPAlignmentMean,randPAlignmentCI,'k','LineWidth',2)
                ax = gca;
                xLim = ax.XLim; yLim = ax.YLim;
                text(0.5,alignT + (yLim(2)-yLim(1))*.1,['p = ',num2str(pT)])
                text(2.5,alignP + (yLim(2)-yLim(1))*.1,['p = ',num2str(pP)])
                ylabel('Alignment Index')
                xticklabels({'P Subspace to T Variability', 'Rand Subspace to T Variability',...
                    'T Subspace to P Variability', 'Rand Subspace to P Variability'});
                xtickangle(45)
                ax.FontSize = 16;
                ax.YLim = [0 1];
                set(gca, 'TickDir', 'out')  
                if saveFig
                    saveas(gcf,fullfile(saveDir,[dataset,'_AlignmentIndex.svg']));
                    saveas(gcf,fullfile(saveDir,[dataset,'_AlignmentIndex.fig']));
                end
        
    end