clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Figure 2';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Turn on/off intermediate figures
    figsOn = true;
    
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
    for datasetList = {'E20200316'}%,'E20200317','E20200318'}%,'N20171215','N20180221','R20201020','R20201021'}%,'R20201020','R20201021'}    %{'E20200317','N20180221'}      
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
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
        end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        
        % Get minimum number of timestamps in condition averages
        numTimestamps = [];
        for i = 1:size(trajStruct,2)
            numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
        end
        if figsOn
            figure
            histogram(numTimestamps)
            xlabel('Number of 25ms bins')
            ylabel('Number of trials')
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
    
            
            
        %Peform LDA by posture on all data 
        obsStruct = struct('label',[],'numObs',[],'allObs',[]);
        structInd = 1;   
        for posture = postureList
           obsStruct(structInd).label = posture;
           allObs = [];
           tempTrajStruct = trajStruct([trajStruct.posture]==posture);
           for i = 1:size(tempTrajStruct,2)
               for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
                   timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
                   traj = tempTrajStruct(i).allSmoothFR(j).traj;
                   if size(traj,1) > minNumTimestamps
                        traj = traj(1:minNumTimestamps,:);
                   end
                   allObs = vertcat(allObs,traj);
               end
           end
           obsStruct(structInd).allObs = allObs;
           obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
           structInd = structInd + 1;
        end
        [postureLDA] = doLDA(obsStruct);
            
            
            
            
            
        if figsOn
            fs = 14;
            %Plot T signal in T subspace    
            figure
            hold on
            for target = targetList
                traj = targetMargTraj((target-1)*minNumTimestamps+1:target*minNumTimestamps,:);
                targetProj = traj*targetPCs;
                plot(targetProj(:,1),targetProj(:,2),'LineWidth',2,'Color','k');
            end
            xlabel('Target Dim 1')
            ylabel('Target Dim 2')
            grid on;
            xticklabels({}); yticklabels({});
            set(gca,'FontSize',fs);
            set(gca,'FontName','Arial');
            %Plot P signal in P subspace
            figure
            hold on
            for posture = postureList
                traj = postureMargTraj((posture-1)*minNumTimestamps+1:posture*minNumTimestamps,:);
                postureProj = traj*posturePCs;
                plot(-postureProj(:,1),postureProj(:,2),'LineWidth',2,'Color',pcmap(posture,:));
            end
            xlabel('Posture Dim 1')
            ylabel('Posture Dim 2')
            grid on;
            xticklabels({}); yticklabels({});
            set(gca,'FontSize',fs);
            set(gca,'FontName','Arial');
        % Visualize data in PT PC subspace
            figure
            hold on
            dims = [1,2,1];
            for posture = postureList
                for target = targetList
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
                    traj = (traj(1:minNumTimestamps,:)-allMu)*allPCs(:,1:numManPCs);
                    targetProj = traj*targetPCs;
                    postureProj = traj*posturePCs;
                    plot3(targetProj(:,dims(1)),targetProj(:,dims(2)),postureProj(:,dims(3)),'LineWidth',2,'Color',pcmap(posture,:));
                end
        %         plot3(targetMeansTProj(:,1),targetMeansTProj(:,2),targetMeansPProj(:,1),'-','LineWidth',2,'Color',pcmap(posture,:));
            end

            figure
            hold on
            dims = [1,2,3];
            for posture = postureList
                for target = targetList
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
                    traj = (traj(1:minNumTimestamps,:)-allMu)*allPCs(:,1:numManPCs);
                    postureProj = traj*posturePCs;
                    plot3(postureProj(:,dims(1)),postureProj(:,dims(2)),postureProj(:,dims(3)),'LineWidth',2,'Color',pcmap(posture,:));
                end
            end

        end
            % Plot cross-projection eigenspectra   
            fs = 20;
            figure
                b = bar([targetTPC_VAF(1:numTargetPCs),postureTPC_VAF(1:numTargetPCs)],'FaceColor','flat');
                for k = 1:2
                   b(k).CData = .7*(k-1)*[1 1 1]; 
                end
                xlabel('Target PC')
                ylabel('Variance Explained (%)')
                legend('Target Signal','Posture Signal')
                set(gca, 'TickDir', 'out')
                set(gca,'fontname','arial')
                set(gca,'fontsize',fs)
                if saveFig
                    saveas(gcf,fullfile(saveDir,[dataset,'_TargetPCSpec.svg']));
                    saveas(gcf,fullfile(saveDir,[dataset,'_TargetPCSpec.fig']));
                end
            figure
                b = bar([targetPPC_VAF(1:numPosturePCs),posturePPC_VAF(1:numPosturePCs)],'FaceColor','flat');
                for k = 1:2
                   b(k).CData = .7*(k-1)*[1 1 1]; 
                end
                xlabel('Posture PC')
                ylabel('Variance Explained (%)')
                set(gca, 'TickDir', 'out')  
                set(gca,'fontname','arial')
                set(gca,'fontsize',fs)
                if saveFig
                    saveas(gcf,fullfile(saveDir,[dataset,'_PosturePCSpec.svg']));
                    saveas(gcf,fullfile(saveDir,[dataset,'_PosturePCSpec.fig']));
                end
        if figsOn
            %Cumulative version of cross-projections
            fs = 20;
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
        end
        
        %How much target-related variance is captured in posture subspace?  What's the most you could hope to capture w an equal number of target dims?
        targetPCaptured = sum(diag(cov(targetMargTraj*posturePCs(:,1:numPosturePCs))));
        targetTCaptured = sum(diag(cov(targetMargTraj*targetPCs(:,1:numPosturePCs))));
    
        %How much posture-related variance is captured in target subspace?  What's the most you could hope to capture w an equal number of posture dims?
        postureTCaptured = sum(diag(cov(postureMargTraj*targetPCs(:,1:numTargetPCs))));
        posturePCaptured = sum(diag(cov(postureMargTraj*posturePCs(:,1:numTargetPCs))));

        %Alignment Index
        Tsig_to_Pspace_AI = targetPCaptured/targetTCaptured;
        Psig_to_Tspace_AI = postureTCaptured/posturePCaptured;
    
        %Generate random subspaces; ask what fraction of target-related and posture-related variability is captured by them
            numDraws = 10000;

            %Target
            randTVarCaptured = zeros(1,numDraws);
            for i = 1:numDraws
               %Draw a random orthonormalized subspace from manifold that has #posturePCs dims
               v = normrnd(0,1,numManPCs,numPosturePCs);
               randomDirs = orth(v);
               %Compute amount of target-related variability that is captured by that space
               randTVarCaptured(i) = sum(diag(cov(targetMargTraj*randomDirs)));
            end
            Tsig_to_Pspace_nullDist = randTVarCaptured/targetTCaptured;

            %Posture
            randPVarCaptured = zeros(1,numDraws);
            for i = 1:numDraws
               %Draw a random orthonormalized subspace from manifold that has #posturePCs dims
               v = normrnd(0,1,numManPCs,numTargetPCs);
               randomDirs = orth(v);
               %Compute amount of posture-related variability that is captured by that space
               randPVarCaptured(i) = sum(diag(cov(postureMargTraj*randomDirs)));
            end
            Psig_to_Tspace_nullDist = randPVarCaptured/posturePCaptured;
    
        %Plot p-value histograms     
        if figsOn
            figure
                hist(Tsig_to_Pspace_nullDist)
                hold on
                ax = gca;
                xLim = ax.XLim; yLim = ax.YLim;
                line([Tsig_to_Pspace_AI,Tsig_to_Pspace_AI],yLim,'Color','r','LineWidth',2);
                xlabel('Alignment Index')
                ylabel('Count')

             figure
                hist(Psig_to_Tspace_nullDist)
                hold on
                ax = gca;
                xLim = ax.XLim; yLim = ax.YLim;
                line([Psig_to_Tspace_AI,Psig_to_Tspace_AI],yLim,'Color','r','LineWidth',2);
                xlabel('Alignment Index')
                ylabel('Count')
        end
        
          %Fill resultStruct
          resultStruct(resultStructInd).animal = dataset(1);
          resultStruct(resultStructInd).dataset = dataset;
          resultStruct(resultStructInd).Psig_to_Tspace_AI = Psig_to_Tspace_AI;
          resultStruct(resultStructInd).Psig_to_Tspace_nullDist = Psig_to_Tspace_nullDist;
          resultStruct(resultStructInd).Tsig_to_Pspace_AI = Tsig_to_Pspace_AI;
          resultStruct(resultStructInd).Tsig_to_Pspace_nullDist = Tsig_to_Pspace_nullDist;
          resultStructInd = resultStructInd + 1;
          
    end

%% Local function for performing LDA - no balancing 
     function [LDAproj] = doLDA(obsStruct)
        numClasses = size(obsStruct,2);
        numDims = size(obsStruct(1).allObs,2);
        obs = []; labels = [];
        for i = 1:numel(obsStruct)
            obs = vertcat(obs,obsStruct(i).allObs);
            labels = vertcat(labels,obsStruct(i).label*ones(size(obsStruct(i).allObs,1),1));
        end
        LDAproj = fisherLDA(obs, labels);
        [LDAproj,~] = qr(LDAproj);
        LDAproj = LDAproj(:,1:numClasses-1);
     end 