clear; clc; clf; close all;
set(0, 'DefaultFigureRenderer', 'painters');

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\Orthogonality Figure';
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;

%% Load Data, subselect
%     task = 'BCI';
%     [Data] = subselectForTrajDist(Data,task);

    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;
    

    %Data = loadData('E20200317');
    %condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
    %trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    %trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
    

    Data = loadData('N20180221');
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
    trialInclStates(1).inclStates = {{'state','Cursor Release','last',-50},{'state','Target Hold','last',0}};
    
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    
%% Get minimum number of timestamps in condition averages
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    figure
    histogram(numTimestamps)
    xlabel('Number of 25ms bins')
    ylabel('Number of trials')
    
    % Get minimum number of trials and timestamps
    [minNumTimestamps,i] = min(numTimestamps);

%% Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);
    
%% Form X containing trial-averaged data for each condition
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
    
%% Do marginalizations of X
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Offsets
    targetMargOffset = mean(X,[1 3],'omitnan') - CIMargOffset;
    postureMargOffset = mean(X,[1 2],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
    targetTrajNoCP = X - CIMargTraj - CIMargOffset - postureMargTraj - postureMargOffset;


%% Do PCA on marginalizations
    allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
    Ca = cov(allTraj);
    [Ua,Sa] = eig(Ca);
    
%     [Ua,score,latent,tsquared,explained,mu] = pca(allTraj); 
    
    CIMargTraj = squeeze(CIMargTraj);
    CIMargOffset = squeeze(CIMargOffset);
    [CIPCA,score,latent,tsquared,explained,mu] = pca(CIMargTraj); 
    
    targetMargTraj = squeeze(targetMargTraj);
    targetMargTraj = reshape(targetMargTraj,[numTargets*minNumTimestamps,numChannels]);
    [Dt,score,latent,tsquared,explained,mu] = pca(targetMargTraj); 
    
    targetTrajNoCP = squeeze(targetTrajNoCP);
    targetTrajNoCP = reshape(targetTrajNoCP,[numTargets*numPostures*minNumTimestamps,numChannels]);
    [targetNoCPPCA,score,latent,tsquared,explained,mu] = pca(targetTrajNoCP); 
      
    postureMargTraj = squeeze(postureMargTraj);
    postureMargTraj = reshape(postureMargTraj,[numPostures*minNumTimestamps,numChannels]);
    [Dp,score,latent,tsquared,explained,mu] = pca(postureMargTraj); 

%% Get cross-projection VAF
    totalTargetVar = trace(cov(targetMargTraj));
    totalPostureVar = trace(cov(postureMargTraj));
    targetTPC_VAF = diag(cov(targetMargTraj*Dt))./totalTargetVar.*100;
    postureTPC_VAF = diag(cov(postureMargTraj*Dt))./totalPostureVar.*100;
    targetPPC_VAF = diag(cov(targetMargTraj*Dp))./totalTargetVar.*100;
    posturePPC_VAF = diag(cov(postureMargTraj*Dp))./totalPostureVar.*100;

%% Keep Enough PCs to capture 90% of respective variance
    varThreshold = 90;
    
    numTargetPCs = 1;
    while sum(targetTPC_VAF(1:numTargetPCs)) < varThreshold
        numTargetPCs = numTargetPCs + 1;
    end
    
    numPosturePCs = 1;
    while sum(posturePPC_VAF(1:numPosturePCs)) < varThreshold
        numPosturePCs = numPosturePCs + 1;
    end
    
%% Visualize data in PT PC subspace
    figure
    hold on
    dims = [1,2,1];
    for posture = postureList
        for target = targetList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
%             traj = traj(1:13,:)-CIMargTraj;
            traj = traj(1:13,:);
            targetProj = traj*Dt;
            postureProj = traj*Dp;
            plot3(targetProj(:,dims(1)),targetProj(:,dims(2)),postureProj(:,dims(3)),'LineWidth',2,'Color',pcmap(posture,:));
        end
%         plot3(targetMeansTProj(:,1),targetMeansTProj(:,2),targetMeansPProj(:,1),'-','LineWidth',2,'Color',pcmap(posture,:));
    end


%% Plot cross-projection eigenspectra   
    figure
        bar([targetTPC_VAF(1:numTargetPCs),postureTPC_VAF(1:numTargetPCs)])
        ax = gca;
        ax.FontSize = 16;
        xlabel('Target PC')
        ylabel('Variance Explained (%)')
        legend('Target Marginalization','Posture Marginalization')
        set(gca, 'TickDir', 'out')
        if saveFig
            saveas(gcf,fullfile(saveDir,'TargetPCSpec.svg'));
            saveas(gcf,fullfile(saveDir,'TargetPCSpec.fig'));
        end
    figure
        bar([targetPPC_VAF(1:numPosturePCs),posturePPC_VAF(1:numPosturePCs)])
        ax = gca;
        ax.FontSize = 16;
        xlabel('Posture PC')
        ylabel('Variance Explained (%)')
        set(gca, 'TickDir', 'out')  
        if saveFig
            saveas(gcf,fullfile(saveDir,'PosturePCSpec.svg'));
            saveas(gcf,fullfile(saveDir,'PosturePCSpec.fig'));
        end
        
%% How much target-related variance is captured in posture subspace?  What's the most you could hope to capture w an equal number of dims?
    targetPCaptured = sum(diag(cov(targetMargTraj*Dp(:,1:numPosturePCs))));
    targetTCaptured = sum(diag(cov(targetMargTraj*Dt(:,1:numPosturePCs))));
    
%% How much posture-related variance is captured in target subspace?  What's the most you could hope to capture w an equal number of dims?
    postureTCaptured = sum(diag(cov(postureMargTraj*Dt(:,1:numTargetPCs))));
    posturePCaptured = sum(diag(cov(postureMargTraj*Dp(:,1:numTargetPCs))));

%% Alignment Index
    alignT = targetPCaptured/targetTCaptured;
    alignP = postureTCaptured/posturePCaptured;
    
%% Generate random subspaces; ask what fraction of target-related and posture-related variability is captured by them
    numDraws = 10000;
    
    %Target
    randTVarCaptured = zeros(1,numDraws);
    for i = 1:numDraws
       %Draw a random orthonormalized subspace that matches data covariance
       %& has #posturePCs dims
       v = normrnd(0,1,87,numPosturePCs);
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
       v = normrnd(0,1,87,numTargetPCs);
       Dr = orth((Ua*sqrt(Sa)*v)/norm(Ua*sqrt(Sa)*v));
       %Compute amount of target-related and posture-related variability that is captured by that
       %space
       randPVarCaptured(i) = sum(diag(cov(postureMargTraj*Dr)));
    end
    randPAlignmentDist = randPVarCaptured/posturePCaptured;
    
%% Plot p-value histograms     
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

%% Plot results
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
        saveas(gcf,fullfile(saveDir,'AlignmentIndex.svg'));
        saveas(gcf,fullfile(saveDir,'AlignmentIndex.fig'));
    end