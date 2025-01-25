clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20221123\Figures\Figure 2\Alignment Indices';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],'pSubspaceAligment',[],'nullMean',[],'nullStd',[]);
    
%% Parameters
    % How many PC's to describe the 'manifold'? 
    manVarThreshold = 80;
    % Keep Enough PCs to capture 90% of respective variance
    varThreshold = 80;
    
%% Run loop for each dataset   
    resultStructInd = 1;
    
%     task = 'BCI';
%     bciDatasetList = {'E20200316','E20200317','E20200318','N20171215','N20180221','R20201020','R20201021'};
%     
%     task = 'Iso';
%     isoDatasetList = {'E20200116','E20200117'};
    
    task = 'Reaching';
    reachDatasetList = {'E20210706','E20210707','E20210708','E20210709','E20210710'...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};
    
    for datasetList = {'E20200316'}%reachDatasetList%{'E20200316','E20200317','E20200318'}%{'E20200316','E20200317','E20200318'}%,'N20171215','N20180221','R20201020','R20201021'}%,'R20201020','R20201021'}    %{'E20200317','N20180221'}      
        %Load data
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
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
            %Iso
            case {'E20200116','E20200117','E20200120'}
                trialInclStates(1).trialName = {'IsometricForce_1D'};   
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',200}};
            %Reaching
            case{'E20210706','E20210707','E20210708','E20210709','E20210710'}
                trialInclStates(1).trialName = {'GridReaching'};  
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
            case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
            case{'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
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
        
        %Keep only postures with every target for Earl
        if strcmpi(dataset(1),'E')
            keepPosture = [];
            for posture = postureList
                tempTrajStruct = trajStruct([trajStruct.posture]==posture);
                postureTargetList = [tempTrajStruct.target];
                if isequal(postureTargetList,targetList)
                    keepPosture = [posture,keepPosture];
                end
            end
            
            trajStruct = trajStruct(ismember([trajStruct.posture],keepPosture));
        end
        
        %Update posture and list
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        
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
    
        %How much target-related variance is captured in posture subspace?  What's the most you could hope to capture w an equal number of target dims?
        targetPCaptured = sum(diag(cov(targetMargTraj*posturePCs(:,1:numPosturePCs))));
        targetTCaptured = sum(diag(cov(targetMargTraj*targetPCs(:,1:numPosturePCs))));
    
        %How much posture-related variance is captured in target subspace?  What's the most you could hope to capture w an equal number of posture dims?
        postureTCaptured = sum(diag(cov(postureMargTraj*targetPCs(:,1:numTargetPCs))));
        posturePCaptured = sum(diag(cov(postureMargTraj*posturePCs(:,1:numTargetPCs))));

        %Alignment Index
        Tsig_to_Pspace_AI = targetPCaptured/targetTCaptured;
        Psig_to_Tspace_AI = postureTCaptured/posturePCaptured;

        %Compute nullDist
        %Generate random subspaces; ask what fraction of target-related and posture-related variability is captured by them
        numDraws = 10000;

        randPVarCaptured = zeros(1,numDraws);
        for i = 1:numDraws
           %Draw a random orthonormalized subspace from manifold that has #posturePCs dims
           v = normrnd(0,1,numManPCs,numTargetPCs);
           randomDirs = orth(v);
           %Compute amount of posture-related variability that is captured by that space
           randPVarCaptured(i) = sum(diag(cov(postureMargTraj*randomDirs)));
        end
        Psig_to_Tspace_nullDist = randPVarCaptured/posturePCaptured;
        
        
      %Fill resultStruct
      resultStruct(resultStructInd).animal = dataset(1);
      resultStruct(resultStructInd).dataset = dataset;
      resultStruct(resultStructInd).pSubspaceAligment = Psig_to_Tspace_AI;
      resultStruct(resultStructInd).nullDist = Psig_to_Tspace_nullDist;
      resultStructInd = resultStructInd + 1;
          
    end


%% Combine within animal; make plot
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).dataset(1);
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
    fs = 16;
    figure; hold on;
    monkeyInd = 1;
    for monkey = monkeyList
       tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
       pSubspaceAligment = horzcat(tempResultStruct.pSubspaceAligment);
       numSamp = size(pSubspaceAligment,2);
       %Plot Data
       scatter(monkeyInd*ones(1,numSamp),pSubspaceAligment,100,'r','filled','MarkerFaceAlpha',0.8);
       nullDist = horzcat(tempResultStruct.nullDist);
       nullMean = mean(nullDist);
       nullStd = std(nullDist);
       %Plot errorbar
       e = errorbar(monkeyInd+.1,nullMean,nullStd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','k');
       alpha = 0.3;
       set([e.Bar, e.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*alpha])
       monkeyInd = monkeyInd + 1;
    end
    
    ax = gca;
    ax.YLim = [0 1]; ax.XLim = [.5 max(monkeyInd-1)+.5];
    yticks([0 1]); xticks(1:max(monkeyInd-1));
    xticklabels(monkeyList)
    set(gca, 'TickDir', 'out')
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    ylabel('PS-TS Alignment')
    xlabel('Monkey')
    if saveFig
        saveas(gcf,fullfile(saveDir,[task,'_alignmentInd.svg']));
    end