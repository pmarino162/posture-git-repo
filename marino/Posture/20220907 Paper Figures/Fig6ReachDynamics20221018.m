clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Figure 8';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    
%% Run loop for each dataset
for datasetList = {'E20210706','N20190228','R20200222'}%{'E20210706'}%,'N20190306','R20200222'}
    %Load data
    dataset = datasetList{1,1};
    [Data,zScoreParams] = loadData(dataset);

    %Get trajStruct
    binWidth = 25; kernelStdDev = 25;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    switch dataset
        case 'E20210706'
            trialInclStates(1).trialName = {'GridReaching'};     
        case 'N20190228'
            trialInclStates(1).trialName = {'Nigel Dissociation'};
        case 'R20200222'
            trialInclStates(1).trialName = {'Rocky Dissociation'};
    end
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    
    %Keep only postures 1-4 for Earl
    if strcmpi(dataset(1),'E')
        trajStruct = trajStruct(ismember([trajStruct.posture],1:4));
    end
    
    % Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);
    
    % Get minimum number of timestamps in condition averages
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    figure
    histogram(numTimestamps)
    xlabel('Number of 25ms bins')
    ylabel('Number of trials')
    [minNumTimestamps,i] = min(numTimestamps);
    
    %Get obsStruct
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    label = 1;
    for posture = postureList
        allObs = [];
        for target = targetList
            allSmoothFR = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allSmoothFR;
            for i = 1:size(allSmoothFR,2)
               obs = allSmoothFR(i).traj;
               allObs = vertcat(obs,allObs);
            end
        end
        obsStruct(label).label = label;
        obsStruct(label).allObs = allObs;
        obsStruct(label).numObs = size(allObs,1);
        label = label + 1;
    end
    
    %Do PCA on observations
    allObs = vertcat(obsStruct.allObs);
    [coeff,score,latent,tsquared,explained,mu] = pca(allObs);
    
    %Replace obsStruct observations w PC projections
    numPCs = 30;
    for i = 1:numel(obsStruct)
       obsStruct(i).allObs = obsStruct(i).allObs*coeff(:,1:numPCs);       
    end
    
    %Do LDA on obsStruct
    [postureLDA] = doLDA(obsStruct);
    
    %Combine PCA and LDA to get posture Axis 
    postureLDA = coeff(:,1:numPCs)*postureLDA;
     
    % Find peri-reach target plane   
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
    allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
    
    %Perform marginalizations and store them
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Offsets
    targetMargOffset = mean(X,[1 3],'omitnan') - CIMargOffset;
    postureMargOffset = mean(X,[1 2],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj - targetMargOffset - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - postureMargOffset - CIMargOffset;
    targetTrajNoCP = X - CIMargTraj - CIMargOffset - postureMargTraj - postureMargOffset;
    postureTrajNoC = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
            
    %Do PCA on Marginalizations    
    CIMargTraj = squeeze(CIMargTraj);
        CIMargOffset = squeeze(CIMargOffset);
        [CIPCA,score,latent,tsquared,explained,mu] = pca(CIMargTraj); 
    targetMargTraj = squeeze(targetMargTraj);
        targetMargTraj = reshape(targetMargTraj,[numTargets*minNumTimestamps,numChannels]);
        [targetPCA,score,latent,tsquared,explained,mu] = pca(targetMargTraj); 
    targetTrajNoCP = squeeze(targetTrajNoCP);
        targetTrajNoCP = reshape(targetTrajNoCP,[numTargets*numPostures*minNumTimestamps,numChannels]);
        [targetNoCPPCA,score,latent,tsquared,explained,mu] = pca(targetTrajNoCP); 
    postureMargOffset = squeeze(postureMargOffset);
        [posturePCA,score,latent,tsquared,explained,mu] = pca(postureMargOffset); 
    postureTrajNoC = squeeze(postureTrajNoC);
        postureTrajNoC = reshape(postureTrajNoC,[numPostures*minNumTimestamps,numChannels]);
        [postureNoCPCA,score,latent,tsquared,explained,mu] = pca(postureTrajNoC); 
      
    %Orthonormalize combinations of axes
%     [CPTOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1),CIPCA(:,1)]); CPTOrth = CPTOrth(:,1:3);
%     [CPOrth,~] = qr([postureLDA(:,1:2),CIPCA(:,1)]); CPOrth = CPOrth(:,1:3);
    [PTOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
        
    %Flip Posture Dim so that posture 1 is on the bottom
    p1TrajMean = mean(trajStruct([trajStruct.posture]==1 & [trajStruct.target]==1).avgSmoothFR.traj,1); 
    p2TrajMean = mean(trajStruct([trajStruct.posture]==max(postureList) & [trajStruct.target]==1).avgSmoothFR.traj,1); 
    p1TrajCoord = p1TrajMean*PTOrth(:,1);
    p2TrajCoord = p2TrajMean*PTOrth(:,1);
    if p1TrajCoord > p2TrajCoord
        PTOrth(:,1) = -1.*PTOrth(:,1);
    end
        
    %Add Projections to trajStruct
    totalVar = trace(cov(allTraj));
    for i = 1:size(trajStruct,2)
        %Add average traces to trajStruct
        trajStruct(i).PTOrth.traj = trajStruct(i).avgSmoothFR.traj*PTOrth;
        trajStruct(i).postureLDA.traj = trajStruct(i).avgSmoothFR.traj*postureLDA;
        %Get VAF
        trajStruct(i).PTOrth.VAF =  100.*(diag(cov(allTraj*PTOrth))')./totalVar;
        trajStruct(i).postureLDA.VAF = 100.*(diag(cov(allTraj*postureLDA))')./totalVar;
    end    
    
    %Make figures
    fs = 14;
    figure; hold on
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(1:minNumTimestamps,:);
                traj = traj-CIMargTraj(1:minNumTimestamps,:);
                traj = traj*PTOrth;
                %traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(1:minNumTimestamps,:);
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end
    grid on
    switch dataset
        case 'E20210706'
            viewVec = [-34.7680,14.0271];
        case 'N20190228'
            viewVec = [50.3536,23.8403];
        case 'R20200222'
            viewVec = [111.7415,34.2224];
    end
    view(viewVec);
    ax = gca;
    xlims = ax.XLim; ylims = ax.YLim; zlims = ax.ZLim;
    xticklabels({}); yticklabels({}); zticklabels({}); 
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    %xlabel('Target Dim 1'); ylabel('Target Dim 2'); zlabel('Posture Dim 1')
    
    VAF = round(trajStruct(1).PTOrth.VAF);
    xlabel(['Target Dim 1 (',num2str(VAF(xDim)),'%)'])
    ylabel(['Target Dim 2 (',num2str(VAF(yDim)),'%)']) 
    zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
    if saveFig
        saveas(gcf,fullfile(saveDir,[dataset,'reachTraj.svg']));
    end
    
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