clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;

%% Load Data, subselect
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;
    task = 'BCI';
    [Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
%     [Data] = subselectForTrajDist(Data,task);
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
    %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-150},{'state','Step 2','first',0}};
    trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
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

%% For each timecourse in trajstruct, get single point
    for i = 1:size(trajStruct,2)
        %Individual trials
        for j = 1:size(trajStruct(i).allSmoothFR,2)
           trajLength = length(trajStruct(i).allSmoothFR(j).timestamps);
           if trajLength > minNumTimestamps
               trajLength = minNumTimestamps;
           end
           trajStruct(i).allSmoothFR(j).trialAvg = mean(trajStruct(i).allSmoothFR(j).traj(1:trajLength,:)); 
        end
        %Condition average
        trajStruct(i).avgSmoothFR.condAvg = mean(trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps,:));
    end

%% Get all avgs; num dims; Ca; Da; etc
    allAvgs = [];
    for i = 1:size(trajStruct,2)
       allAvgs = vertcat(allAvgs,trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps,:)); 
    end
    numDims = size(allAvgs,2);
    Ca = cov(allAvgs);
    [Ua,Sa] = eig(Ca);
    [Da,~,~,~,Expa,mua] = pca(allAvgs);
    
%% Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    
%% Do PCA on all means
    allMeans = zeros(40,87);
    for i = 1:size(trajStruct,2)
       allMeans(i,:) = trajStruct(i).avgSmoothFR.condAvg; 
    end
    [DallMeans,~,~,~,ExpAllMeans,muAllMeans] = pca(allMeans);
    
    figure
    hold on
    dims = [1,2,3];
    for posture = postureList
        targetMeans = zeros(numTargets,87);
        for target = targetList
            targetMeans(target,:) = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.condAvg;
        end
        targetMeansProj = targetMeans*DallMeans;
        for target = targetList
            plot3(targetMeansProj(target,dims(1)),targetMeansProj(target,dims(2)),targetMeansProj(target,dims(3)),'.','MarkerSize',20,'Color',tcmap(target,:));
        end
        plot3(targetMeansProj(:,dims(1)),targetMeansProj(:,dims(2)),targetMeansProj(:,dims(3)),'-','LineWidth',2,'Color',pcmap(posture,:));
    end
    axis equal
    
%% Get Posture Means; Posture subspace
    postureMeans = zeros(numPostures,numDims,numTargets);
    i = 1;
    for posture = postureList
        tempData = trajStruct([trajStruct.posture]==posture);
        for j = 1:size(tempData,2)
            postureMeans(i,:,j) = tempData(j).avgSmoothFR.condAvg;
        end
        i = i+1;
    end
    postureMeans = mean(postureMeans,3);
    [Dp,~,~,~,Expp,mup] = pca(postureMeans); 
    
%% Get Target Means; Target Subspace
    targetMeans = zeros(numPostures,numDims,numTargets);
    i = 1;
    for target = targetList
        tempData = trajStruct([trajStruct.target]==target);
        for j = 1:size(tempData,2)
            targetMeans(i,:,j) = tempData(j).avgSmoothFR.condAvg;
        end
        i = i+1;
    end
    targetMeans = mean(targetMeans,3);
    [Dt,~,~,~,Expt,mut] = pca(targetMeans); 
       
%% Visualize data in PT PC subspace 
    figure
    hold on
    dims = [1,2,3];
    for posture = postureList
        targetMeans = zeros(numTargets,87);
        for target = targetList
            targetMeans(target,:) = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.condAvg;
        end
        targetMeansTProj = targetMeans*Dt;
        targetMeansPProj = targetMeans*Dp;
        for target = targetList
            plot3(targetMeansTProj(target,1),targetMeansTProj(target,2),targetMeansPProj(target,1),'.','MarkerSize',20,'Color',tcmap(target,:));
        end
        plot3(targetMeansTProj(:,1),targetMeansTProj(:,2),targetMeansPProj(:,1),'-','LineWidth',2,'Color',pcmap(posture,:));
    end

%% Get cross-projection VAF
    totalTargetVar = trace(cov(targetMeans));
    totalPostureVar = trace(cov(postureMeans));
    targetTSubspaceVAF = diag(cov(targetMeans*Dt))./totalTargetVar.*100;
    postureTSubspaceVAF = diag(cov(postureMeans*Dt))./totalPostureVar.*100;
    targetPSubspaceVAF = diag(cov(targetMeans*Dp))./totalTargetVar.*100;
    posturePSubspaceVAF = diag(cov(postureMeans*Dp))./totalPostureVar.*100;

    figure
        bar([targetTSubspaceVAF,postureTSubspaceVAF])
        xlabel('Target PC')
        ylabel('Variance Explained (%)')
        legend('Target Means','Posture Means')

    figure
        bar([targetPSubspaceVAF,posturePSubspaceVAF])
        xlabel('Posture PC')
        ylabel('Variance Explained (%)')
        legend('Target Means','Posture Means')
    
%% Keep first 2 dims of target and posture subspaces
    Dt = Dt(:,1:2);
    Dp = Dp(:,1:2);
    
%% How much target-related variance is captured by target and posture subspaces?
    targetTCaptured = sum(diag(cov(targetMeans*Dt)));
    targetPCaptured = sum(diag(cov(targetMeans*Dp)));
    
%% How much posture-related variance is captured by target and posture subspaces?
    postureTCaptured = sum(diag(cov(postureMeans*Dt)));
    posturePCaptured = sum(diag(cov(postureMeans*Dp)));

%% Alignment Index
    alignT = targetPCaptured/targetTCaptured;
    alignP = postureTCaptured/posturePCaptured;
    
%% Generate random 2-d subspaces; ask what fraction of target-related and posture-related variability is captured by them
    numDraws = 10000;
    randTVarCaptured = zeros(1,numDraws);
    randPVarCaptured = zeros(1,numDraws);
    for i = 1:numDraws
       %Draw a two-D orthonormalized subspace that matches data covariance
       v = normrnd(0,1,87,2);
       Dr = orth((Ua*sqrt(Sa)*v)/norm(Ua*sqrt(Sa)*v));
       %Compute amount of target-related and posture-related variability that is captured by that
       %space
       randTVarCaptured(i) = sum(diag(cov(targetMeans*Dr)));
       randPVarCaptured(i) = sum(diag(cov(postureMeans*Dr)));
    end
    randTAlignmentDist = randTVarCaptured/targetTCaptured;
    randPAlignmentDist = randPVarCaptured/posturePCaptured;
    
    figure
        hist(randTAlignmentDist)
        hold on
        ax = gca;
        xLim = ax.XLim; yLim = ax.YLim;
        line([alignT,alignT],yLim,'Color','r','LineWidth',2);
        p = sum(randTAlignmentDist <= alignT)/numDraws;
        text(mean(xLim),yLim(2)-(yLim(2)-yLim(1))*.1,['p = ',num2str(p)])
        xlabel('Alignment Index')
        ylabel('Count')
        
     figure
        hist(randPAlignmentDist)
        hold on
        ax = gca;
        xLim = ax.XLim; yLim = ax.YLim;
        line([alignP,alignP],yLim,'Color','r','LineWidth',2);
        p = sum(randPAlignmentDist <= alignP)/numDraws;
        text(mean(xLim),yLim(2)-(yLim(2)-yLim(1))*.1,['p = ',num2str(p)])
        xlabel('Alignment Index')
        ylabel('Count')
        
    
%     randTVarCapturedMean = mean(randTVarCaptured);
%     randTVarCapturedCI = 1.96*std(randTVarCaptured)/sqrt(numDraws);

    
%% Get RSV for subspaces; compare to random vectors drawn from within data cov
    totalVar = trace(cov(allAvgs));
    Ca = cov(allAvgs);
    
    RSVt = trace(Ca*Dt*Dt')/totalVar*100;
    RSVp = trace(Ca*Dp*Dp')/totalVar*100;
    RSVpt = trace(Ca*Dp*Dp'*Dt*Dt')/totalVar*100;
    RSVtp = trace(Ca*Dt*Dt'*Dp*Dp')/totalVar*100;
    
    numDraws = 10000;
    RSVrand = zeros(1,numDraws);
    for i = 1:numDraws
       %Draw two, two-D orthonormalized subspaces that match data
       %covariance
       v1 = normrnd(0,1,87,2);
       v2 = normrnd(0,1,87,2);
       space1 = orth((Ua*sqrt(Sa)*v1)/norm(Ua*sqrt(Sa)*v1));
       space2 = orth((Ua*sqrt(Sa)*v2)/norm(Ua*sqrt(Sa)*v2));
       %Compute RSV for those spaces, store
       RSVrand(i) = trace(Ca*space1*space1'*space2*space2')/totalVar*100;
    end
    RSVrandmean = mean(RSVrand);
    RSVCI = 1.96*std(RSVrand)/sqrt(numDraws);
    
%% Do statistics 
    figure
        hist(RSVrand)
        hold on
        ax = gca;
        xLim = ax.XLim; yLim = ax.YLim;
        line([RSVtp,RSVtp],yLim,'Color','r','LineWidth',2);
        p = sum(RSVrand <= RSVtp)/numDraws;
        text(mean(xLim),yLim(2)-(yLim(2)-yLim(1))*.1,['p = ',num2str(p)])
        xlabel('RSV (%)')
        ylabel('Count')
    
%% Plot results
    f = figure;
    f.Position = [100 100 500 500];
    bar([RSVtp,RSVrandmean])
    hold on
    errorbar(2,RSVrandmean,RSVCI,'k','LineWidth',2)
    ax = gca;
    xLim = ax.XLim; yLim = ax.YLim;
    text(1,RSVtp + (yLim(2)-yLim(1))*.1,['p = ',num2str(p)])
    ylabel('Retained Subspace Variance (%)')
    xticklabels({'Posture and Target','Random'})
    xtickangle(45)
    ax.FontSize = 16;