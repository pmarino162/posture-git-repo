clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Postdoc applications\Job Talk\Comp Model Animations';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%% Construct resultStruct for multiple sessions
    %Earl BCI
    inputStruct(1).dataset = 'E20200317'; inputStruct(1).task = 'BCI'; inputStruct(1).epoch = 'all';

    resultStruct = [];
    for i = 1:numel(inputStruct)
        dataset = inputStruct(i).dataset;
        task = inputStruct(i).task;
        epoch = inputStruct(i).epoch;
        [tempResultStruct] = leaveOneConditionOutCompAnalysis(dataset,task,epoch);
        resultStruct = vertcat(resultStruct,tempResultStruct);
    end
    
%% Visualize results for a dataset
    %Get actual and predicted traj
    tempResultStruct = resultStruct(1);
    actualTrajStruct = tempResultStruct.actualTrajStruct;
    predTrajStruct = tempResultStruct.predTrajStruct;
    
    % Get minimum number of timestamps in condition averages
    minNumTimestamps = size(actualTrajStruct(1).traj,1);  
    
    %Get posture and target lists
    postureList = unique([actualTrajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([actualTrajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(actualTrajStruct(1).traj,2);
    numConditions = size(actualTrajStruct,2);
    
    %Collect both actual and predicted data into X for PCA and marginalization 
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    X = NaN(minNumTimestamps,2*numTargets,numPostures,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).traj;
            predTraj = predTrajStruct([predTrajStruct.posture]==posture & [predTrajStruct.target]==target).traj;
            X(:,targetInd,postureInd,:) = actualTraj(1:minNumTimestamps,:); 
            X(:,targetInd+1,postureInd,:) = predTraj(1:minNumTimestamps,:); 
            targetInd = targetInd + 2;
        end
        postureInd = postureInd+1;
    end
    
    %Get PC space
    [allPCA,score,latent,tsquared,explained,allMu] = pca(reshape(X,[2*numTargets*numPostures*minNumTimestamps,numChannels]));
    
    %Get CI Space
    CIMargTraj = mean(X,[2 3],'omitnan');
    [CIPCA,score,latent,tsquared,explained,mu] = pca(squeeze(CIMargTraj)); 
    
    %Get target space
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj;
    targetMargTraj = squeeze(targetMargTraj);
    targetMargTraj = reshape(targetMargTraj,[2*numTargets*minNumTimestamps,numChannels]);
    [targetPCA,score,latent,tsquared,explained,mu] = pca(targetMargTraj); 
    
    %Get posture space
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = postureList
       obsStruct(structInd).label = posture;
       allObs = [];
       tempTrajStruct = actualTrajStruct([actualTrajStruct.posture]==posture);
       for i = 1:size(tempTrajStruct,2)
           traj = tempTrajStruct(i).traj;
           allObs = vertcat(allObs,traj);
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    [postureLDA] = doLDA(obsStruct);
    
    %Orthonormalize combinations of axes
    [PTOrth,~] = qr([postureLDA(:,1),targetPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);

    %Flip Posture Dim so that posture 1 is on the left
    p1TrajMean = mean(actualTrajStruct([actualTrajStruct.posture]==1 & [actualTrajStruct.target]==1).traj,1); 
    p2TrajMean = mean(actualTrajStruct([actualTrajStruct.posture]==max(postureList) & [actualTrajStruct.target]==1).traj,1); 
    p1TrajCoord = p1TrajMean*postureLDA(:,1);
    p2TrajCoord = p2TrajMean*postureLDA(:,1);
    if p1TrajCoord > p2TrajCoord
        postureLDA(:,1) = -1.*postureLDA(:,1);
    end
    p1TrajCoord = p1TrajMean*PTOrth(:,1);
    p2TrajCoord = p2TrajMean*PTOrth(:,1);
    if p1TrajCoord > p2TrajCoord
        PTOrth(:,1) = -1.*PTOrth(:,1);
    end

    %Add Projections to trajStruct
    %totalVar = trace(cov(allTraj));
    for i = 1:size(actualTrajStruct,2)
        %Add average traces to trajStruct
        actualTrajStruct(i).PCA = (actualTrajStruct(i).traj)*allPCA;
        predTrajStruct(i).PCA = (predTrajStruct(i).traj)*allPCA;       
        actualTrajStruct(i).PTOrth = (actualTrajStruct(i).traj)*PTOrth;
        predTrajStruct(i).PTOrth = (predTrajStruct(i).traj)*PTOrth; 
    end
        
%% Fig 1 - actual data  Get axis limits.
    fs = 14;
    timePts = 1:12;
    figure; hold on
    tempTrajStruct = actualTrajStruct;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
                traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth(timePts,:); 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end    
    grid on
    xticklabels({}); yticklabels({}); zticklabels({});    
    %xlabel('Target Dim 1'); ylabel('Target Dim 2') ;zlabel('Posture Dim 1')
    view([42.9481 33.9444])
    ax = gca; 
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
   if saveFig
      saveas(gcf,fullfile(saveDir,[dataset,'_ActualPTOrth.svg']));
   end
       xlimits = ax.XLim; ylimits = ax.YLim; zlimits = ax.ZLim;
    xtickVals = ax.XTick; ytickVals = ax.YTick; ztickVals = ax.ZTick;   
    
    
%% Fig 2 - same as fig 1, but turn held out traj purple
purp = tcmap(8,:);
    fs = 14;
    timePts = 1:12;
    figure; hold on
    tempTrajStruct = actualTrajStruct;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth(timePts,:); 
            if posture == 1 && target==4
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',purp,'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',purp,'MarkerFaceColor',purp);
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',purp,'MarkerFaceColor',purp);
            else
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end    
    grid on
    xticklabels({}); yticklabels({}); zticklabels({});    
    xlim(xlimits); ylim(ylimits); zlim(zlimits);
    xticks(xtickVals); yticks(ytickVals); zticks(ztickVals);
    view([42.9481 33.9444])
    ax = gca; 
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
   if saveFig
      saveas(gcf,fullfile(saveDir,[dataset,'_HeldOutTrajPurple.svg']));
   end   
   
%% Fig 3 - turn other T traj purple; only display them and posture 1
    fs = 14;
    timePts = 1:12;
    figure; hold on
    tempTrajStruct = actualTrajStruct;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth(timePts,:); 
            if posture > 1 && target==4
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',purp,'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',purp,'MarkerFaceColor',purp);
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',purp,'MarkerFaceColor',purp);
            elseif posture == 1 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end    
    grid on
    xticklabels({}); yticklabels({}); zticklabels({});    
%xlabel('Target Dim 1'); ylabel('Target Dim 2') ;zlabel('Posture Dim 1')
    view([42.9481,33.9444])
    ax = gca; 
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
   if saveFig
      saveas(gcf,fullfile(saveDir,[dataset,'_PurpleTTraj.svg']));
   end
 
   
%% Fig 4 - turn IC's purple; only show posture 1
 fs = 14;
    timePts = 1:12;
    figure; hold on
    tempTrajStruct = actualTrajStruct;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth(timePts,:); 
            if posture == 1 && target~=4
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',purp,'MarkerFaceColor',purp);
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            else
            end
        end
    end    
    grid on
    xticklabels({}); yticklabels({}); zticklabels({});    
    %xlabel('Target Dim 1'); ylabel('Target Dim 2') ;zlabel('Posture Dim 1')
    view([42.9481,33.9444])
        xlim(xlimits); ylim(ylimits); zlim(zlimits);
    xticks(xtickVals); yticks(ytickVals); zticks(ztickVals);
    ax = gca; 
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
   if saveFig
      saveas(gcf,fullfile(saveDir,[dataset,'_PurpleICs.svg']));
   end
   
   
%% Fig 1 - actual data, before rotation.  Get axis limits.
    fs = 14;
    timePts = 1:12;
    figure; hold on
    tempTrajStruct = actualTrajStruct;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
                traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth(timePts,:); 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end    
    grid on
    xticklabels({}); yticklabels({}); zticklabels({});    
    %xlabel('Target Dim 1'); ylabel('Target Dim 2') ;zlabel('Posture Dim 1')
    view([30 20])
    ax = gca; 
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
   if saveFig
      saveas(gcf,fullfile(saveDir,[dataset,'_ActualPTOrth.svg']));
   end
   
%% Fig 2 - rotate down
    fileName = fullfile(saveDir,'RotatingDownToShowTraj.mp4');
%     f = figure; set(gcf,'color','white'); hold on
%     ax = gca; ax.XLim = postureXLimits; ax.YLim = postureXLimits; ax.ZLim = zlimits;
    OptionZ.FrameRate=12;OptionZ.Duration=5.5;OptionZ.Periodic=false;
    viewVec = [30,20; 42.9481,33.9444];
    CaptureFigVid(viewVec,fileName,OptionZ)
    
%% Fig 3 - turn other T traj purple
purp = tcmap(8,:);

    fs = 14;
    timePts = 1:12;
    figure; hold on
    tempTrajStruct = actualTrajStruct;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth(timePts,:); 
            if posture > 1 && target==4
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',purp,'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',purp,'MarkerFaceColor',purp);
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',purp,'MarkerFaceColor',purp);
            else
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end    
    grid on
    xticklabels({}); yticklabels({}); zticklabels({});    
%xlabel('Target Dim 1'); ylabel('Target Dim 2') ;zlabel('Posture Dim 1')
    view([42.9481,33.9444])
    ax = gca; 
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
   if saveFig
      saveas(gcf,fullfile(saveDir,[dataset,'_PurpleTTraj.svg']));
   end

    xlimits = ax.XLim; ylimits = ax.YLim; zlimits = ax.ZLim;
    xtickVals = ax.XTick; ytickVals = ax.YTick; ztickVals = ax.ZTick;
    
%% Fig 4 - display only T avg
        fs = 14;
    timePts = 1:12;
    figure; hold on
    tempTrajStruct = predTrajStruct;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = 1
        for target = 4
            traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth(timePts,:); 
            traj = traj + [1.5 0 0];
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',purp,'LineWidth',2);
            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',purp,'MarkerFaceColor',purp);
            plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',purp,'MarkerFaceColor',purp);
        end
    end    

    xticklabels({}); yticklabels({}); zticklabels({});    
%xlabel('Target Dim 1'); ylabel('Target Dim 2') ;zlabel('Posture Dim 1')
    view([42.9481,33.9444])
    xlim(xlimits); ylim(ylimits); zlim(zlimits);
    xticks(xtickVals); yticks(ytickVals); zticks(ztickVals);
    grid on
    ax = gca; 
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
   if saveFig
      saveas(gcf,fullfile(saveDir,[dataset,'_OnlyTTraj.svg']));
   end



%% Fig 5 - turn IC's purple
 fs = 14;
    timePts = 1:12;
    figure; hold on
    tempTrajStruct = actualTrajStruct;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth(timePts,:); 
            if posture == 1 && target~=4
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',purp,'MarkerFaceColor',purp);
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            else
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end    
    grid on
    xticklabels({}); yticklabels({}); zticklabels({});    
%xlabel('Target Dim 1'); ylabel('Target Dim 2') ;zlabel('Posture Dim 1')
    view([42.9481,33.9444])
    ax = gca; 
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
   if saveFig
      saveas(gcf,fullfile(saveDir,[dataset,'_PurpleICs.svg']));
   end
   
%% Fig 6 - plot model prediction
        fs = 14;
    timePts = 1:12;
    figure; hold on
    tempTrajStruct = predTrajStruct;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = 1
        for target = 4
            traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth(timePts,:); 
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',purp,'LineWidth',2);
            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',purp,'MarkerFaceColor',purp);
            plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',purp,'MarkerFaceColor',purp);
        end
    end    

    xticklabels({}); yticklabels({}); zticklabels({});    
%xlabel('Target Dim 1'); ylabel('Target Dim 2') ;zlabel('Posture Dim 1')
    view([42.9481,33.9444])
    xlim(xlimits); ylim(ylimits); zlim(zlimits);
    xticks(xtickVals); yticks(ytickVals); zticks(ztickVals);
    grid on
    ax = gca; 
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
   if saveFig
      saveas(gcf,fullfile(saveDir,[dataset,'_ModelPrediction.svg']));
   end




%% Fig 7 - rotate back
    fs = 14;
    timePts = 1:12;
    figure; hold on
    tempTrajStruct = actualTrajStruct;
    xDim = 2; yDim = 3; zDim = 1;
    for posture = postureList
        for target = targetList
            if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
                traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth(timePts,:); 
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
    end    
    grid on
    xticklabels({}); yticklabels({}); zticklabels({});    

    view([42.9481,33.9444])
    ax = gca; 
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)

       fileName = fullfile(saveDir,['RotatingBackUp.mp4']);
%     f = figure; set(gcf,'color','white'); hold on
%     ax = gca; ax.XLim = postureXLimits; ax.YLim = postureXLimits; ax.ZLim = zlimits;
    OptionZ.FrameRate=12;OptionZ.Duration=5.5;OptionZ.Periodic=false;
    viewVec = [ 42.9481,33.9444; 30,20;];
    CaptureFigVid(viewVec,fileName,OptionZ)
%% Local function for performing LDA
     function [LDAproj] = doLDA(obsStruct)
        %Preallocate
        minNumObs = min([obsStruct.numObs]);
        numClasses = size(obsStruct,2);
        numDims = size(obsStruct(1).allObs,2);
        obs = NaN(minNumObs*numClasses,numDims); labels =  NaN(minNumObs*numClasses,1); 

        %Fill
        k = 1;
        for i = 1:size(obsStruct,2)
            totalClassObs = size(obsStruct(i).allObs,1);
            obs(k:k+minNumObs-1,:) = obsStruct(i).allObs(randsample(totalClassObs,minNumObs),:);
            labels(k:k+minNumObs-1,1) = ones(minNumObs,1).*obsStruct(i).label;
            k = k+minNumObs;
        end

        LDAproj = fisherLDA(obs, labels);
        [LDAproj,~] = qr(LDAproj);
        LDAproj = LDAproj(:,1:numClasses-1);
     end  