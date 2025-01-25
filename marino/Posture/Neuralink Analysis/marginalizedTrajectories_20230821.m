clear; clc; clf; close all

%% Parameters
    numManPCs = 10;
    
%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Job Applications\Neuralink\Interview Presentation\Separable effects of posture and goal';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\purpleAndGreen.mat')
    tcmap = purpleAndGreen;
    
%% Load data, get trajStruct
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
            trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);


  %Get minimum number of timestamps in condition averages
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    [minNumTimestamps,i] = min(numTimestamps);
    numPts = 9;
    
     %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);


    %Project all data down to top PCs
    allTraj = NaN(numConditions*numPts,numChannels);
    j = 1;
    for i = 1:numConditions
       allTraj(j:j+numPts-1,:) = trajStruct(i).avgSmoothFR.traj(1:numPts,:);
       j = j + numPts;
    end
    [allPCs,~,~,~,explained,allMu] = pca(allTraj);
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgSmoothFR.traj = (trajStruct(i).avgSmoothFR.traj-allMu)*allPCs(:,1:numManPCs); 
       for j = 1:size(trajStruct(i).allSmoothFR,2)
            trajStruct(i).allSmoothFR(j).traj = (trajStruct(i).allSmoothFR(j).traj-allMu)*allPCs(:,1:numManPCs); 
       end
    end
    
    %Form X
    X = NaN(numPts,numTargets,numPostures,numManPCs);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            X(:,targetInd,postureInd,:) = traj(1:numPts,:); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end

%% Perform marginalizations and store them
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
    
    %Interaction
    intMargOffset = mean(X,[1],'omitnan') - targetMargOffset - postureMargOffset - CIMargOffset ;
    intMargTraj = X-intMargOffset-targetMargTraj-postureMargTraj-targetMargOffset-postureMargOffset-CIMargTraj-CIMargOffset;

%% Fit Marg PCs for eigenspectra    
     targetMargTrajReshape = squeeze(targetMargTraj);
        targetMargTrajReshape = reshape(targetMargTrajReshape,[numTargets*numPts,numManPCs]);
        [targetPCs,~,~,~,targetExplained,targetMargMu] = pca(targetMargTrajReshape); 

        figure 
        bar(targetExplained(1:numManPCs,1))
        xlabel('PC')
        ylabel('Variance Explained (%)')
        ax = gca;
        ax.FontName = 'arial';
        ax.FontSize = 14;
        if saveFig
            saveas(gcf,fullfile(saveDir,'targetEigenSpec.svg'));
        end
        
    postureMargTrajReshape = squeeze(postureMargTraj);
        postureMargTrajReshape = reshape(postureMargTrajReshape,[numPostures*numPts,numManPCs]);
        [posturePCs,~,~,~,postureExplained,postureMargMu] = pca(postureMargTrajReshape); 

        figure 
        bar(postureExplained(1:numManPCs,1))
        xlabel('PC')
        ylabel('Variance Explained (%)')
        ax = gca;
        ax.FontName = 'arial';
        ax.FontSize = 14;
        if saveFig
            saveas(gcf,fullfile(saveDir,'postureEigenSpec.svg'));
        end

%% Get VAFs of marginalizations
    %Get total var
    totalVar = 0;
    numObs = 0;
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:numPts
                totalVar = totalVar + (X(i,targetInd,postureInd,:)-CIMargOffset).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    totalVar = sum(totalVar)./(numObs-1);
    
    %Get CI Var 
    CIVar = 0;
    numObs = 0;
    tempMean = mean(CIMargTraj,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:numPts
%                 CIVar = CIVar + (CIMargTraj(i,1,1,:)+CIMargOffset-CIMargOffset).^2;
                CIVar = CIVar + (CIMargTraj(i,1,1,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    CIVar = sum(CIVar)./(numObs-1);
    CIVarPct = (CIVar/totalVar)*100;
    
    %Get Target Var
    targetVar = 0;
    numObs = 0;
    tempMean = mean(targetMargTraj+targetMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:numPts
%                 targetVar = targetVar + (targetMargTraj(i,target,1,:)+targetMargOffset(1,target,1,:)+CIMargTraj(i,1,1,:)+CIMargOffset-CIMargOffset).^2;
                targetVar = targetVar + (targetMargTraj(i,targetInd,1,:)+targetMargOffset(1,targetInd,1,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    targetVar = sum(targetVar)./(numObs-1);
    targetVarPct = (targetVar/totalVar)*100;
    
    
    %Get Target Offset Var
    targetOffsetVar = 0;
    numObs = 0;
    tempMean = mean(targetMargTraj+targetMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:numPts
                targetOffsetVar = targetOffsetVar + (targetMargOffset(1,targetInd,1,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    targetOffsetVar = sum(targetOffsetVar)./(numObs-1);
    targetOffsetVarPct = (targetOffsetVar/targetVar)*100;
    
    
    %Get Posture Var
    postureVar = 0;
    numObs = 0;
    tempMean = mean(postureMargTraj+postureMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:numPts
                %postureVar = postureVar + (postureMargTraj(i,1,posture,:)+postureMargOffset(1,1,posture,:)+CIMargTraj(i,1,1,:)+CIMargOffset-CIMargOffset).^2;
                postureVar = postureVar + (postureMargTraj(i,1,postureInd,:)+postureMargOffset(1,1,postureInd,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    postureVar = sum(postureVar)./(numObs-1);
    postureVarPct = (postureVar/totalVar)*100;
    
    %Get Posture Offset Var
    postureOffsetVar = 0;
    numObs = 0;
    tempMean = mean(postureMargTraj+postureMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:numPts
                postureOffsetVar = postureOffsetVar + (postureMargOffset(1,1,postureInd,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    postureOffsetVar = sum(postureOffsetVar)./(numObs-1);
    postureOffsetVarPct = (postureOffsetVar/postureVar)*100;
    
    %Get Interaction Var
    intVar = 0;
    numObs = 0;
    tempMean = mean(intMargTraj+intMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:numPts
%                 intVar = intVar + (intMargTraj(i,target,posture,:)+intMargOffset(1,target,posture,:)+targetMargTraj(i,target,1,:)+targetMargOffset(1,target,1,:)+postureMargTraj(i,1,posture,:)+postureMargOffset(1,1,posture,:)+CIMargTraj(i,1,1,:)+CIMargOffset-CIMargOffset).^2;
                intVar = intVar + (intMargTraj(i,targetInd,postureInd,:)+intMargOffset(1,targetInd,postureInd,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    intVar = sum(intVar)./(numObs-1);
    intVarPct = (intVar/totalVar)*100;
    
    %Get Int Offset Var Pct
    intOffsetVar = 0;
    numObs = 0;
    tempMean = mean(intMargTraj+intMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:numPts
%                 intVar = intVar + (intMargTraj(i,target,posture,:)+intMargOffset(1,target,posture,:)+targetMargTraj(i,target,1,:)+targetMargOffset(1,target,1,:)+postureMargTraj(i,1,posture,:)+postureMargOffset(1,1,posture,:)+CIMargTraj(i,1,1,:)+CIMargOffset-CIMargOffset).^2;
%                 intVar = intVar + (intMargTraj(i,targetInd,postureInd,:)+intMargOffset(1,targetInd,postureInd,:)-tempMean).^2;
                intOffsetVar = intOffsetVar + (intMargOffset(1,targetInd,postureInd,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    intOffsetVar = sum(intOffsetVar)./(numObs-1);
    intOffsetVarPct = (intOffsetVar/intVar)*100;
    
%% Visualize by projecting marginalizations into PC's
    
    time = trajStruct(1).avgSmoothFR.timestamps(1:numPts);
    
    %Condition-invariant
    fs = 14;
    f = figure; f.Position = [50 50 1000 400];
    traj = squeeze(CIMargTraj);
    for i = 1:6
        subplot(2,3,i); hold on;
        plot(time,traj(:,i),'LineWidth',2)
        ax = gca;
        ax.YLim = [-5 5];
        ylabel(['PC ',num2str(i)])
    end
    for dim = 1:6
        subplot(2,3,dim);
        ax = gca; 
        yticks([min(ax.YTick) max(ax.YTick)]);
        if dim < 4
            xticks([])
        else
            xticks([min(ax.XTick) max(ax.XTick)]);
        end
        if dim == 5
            xlabel('time (ms)')
        end
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
    end      
    sgtitle('Condition-Invariant Component')
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'TimeMarginalization.jpg')); 
    end

    
    %Target
    f = figure; f.Position = [50 50 1000 400];
    allTraj = targetMargTraj+targetMargOffset;
    targetInd = 1;
    for target = targetList
       traj = squeeze(allTraj(:,targetInd,:,:));
       for i = 1:6
            subplot(2,3,i)
            plot(time,traj(:,i),'Color',tcmap(target,:),'LineWidth',2)
            ax = gca;
            ax.YLim = [-5 5];
            hold on;
            ylabel(['PC ',num2str(i)])
       end
       targetInd = targetInd + 1;
    end
    for dim = 1:6
        subplot(2,3,dim);
        ax = gca; 
        yticks([min(ax.YTick) max(ax.YTick)]);
        if dim < 4
            xticks([])
        else
            xticks([min(ax.XTick) max(ax.XTick)]);
        end
        if dim == 5
            xlabel('time (ms)')
        end
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
    end      
    sgtitle('Goal Component')
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'TargetMarginalization.jpg')); 
    end
    
    %Posture
    f = figure; f.Position = [50 50 1000 400];
    allTraj = postureMargTraj+postureMargOffset;
    postureInd = 1;
    for posture = postureList
       traj = squeeze(allTraj(:,:,postureInd,:));
       for i = 1:6
            subplot(2,3,i)
            plot(time,traj(:,i),'Color',pcmap(posture,:),'LineWidth',2)
            ax = gca;
            ax.YLim = [-5 5];
            hold on;
            ylabel(['PC ',num2str(i)])
       end
       postureInd = postureInd + 1;
    end
    for dim = 1:6
        subplot(2,3,dim);
        ax = gca; 
        yticks([min(ax.YTick) max(ax.YTick)]);
        if dim < 4
            xticks([])
        else
            xticks([min(ax.XTick) max(ax.XTick)]);
        end
        if dim == 5
            xlabel('time (ms)')
        end
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
    end      
    sgtitle('Posture Component')
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'PostureMarginalization.jpg')); 
    end
    
  
    