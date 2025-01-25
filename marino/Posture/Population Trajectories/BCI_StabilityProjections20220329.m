clear; clc; clf; close all;
set(0, 'DefaultFigureRenderer', 'painters');
 
%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\MultiTask Projections figure';

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;
    
%% Run loop for each task
    taskCell ={'BCI'}
    task = taskCell{1,1};
    
    %Load Data, subselect, get traj struct
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;
    switch task
        case 'BCI'    
            [Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
            %[Data] = subselectForTrajDist(Data,task);
            trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
            condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%             trajStruct = trajStruct(ismember([trajStruct.target],[1,3,5,7]));
    end
    

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
    
    %Get posture and target lists
    switch task
        case {'BCI','planning','iso','multijoint BCI'}
            postureList = unique([trajStruct.posture]);
        case 'reaching'
            postureList = unique([trajStruct.posture]);
            %postureList = unique([fourPostureTrajStruct.posture]);
            %postureList = unique([planningTrajStruct.posture]);
    end
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);
    
    %Get pcmap for number of postures
    pcmap = parula(numPostures);
    lightpcmap = rgb2hsv(pcmap);
    lightpcmap(:,2)=.3;
    lightpcmap =hsv2rgb(lightpcmap);
    
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
    
    %Peform LDA by posture on all data 
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = postureList
       obsStruct(structInd).label = posture;
       allObs = [];
       switch task
           case 'BCI'
               tempTrajStruct = trajStruct([trajStruct.posture]==posture);
       end
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

    %Orthonormalize combinations of axes
    [CPTOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1),CIPCA(:,1)]); CPTOrth = CPTOrth(:,1:3);
    [CPOrth,~] = qr([postureLDA(:,1:2),CIPCA(:,1)]); CPOrth = CPOrth(:,1:3);
    [PTOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
    [CPTAllOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1:2),CIPCA(:,1)]); CPTAllOrth = CPTAllOrth(:,1:4);
    
    %Flip C dim if necessary
    traj = trajStruct(1).avgSmoothFR.traj;
    CPT = trajStruct(1).avgSmoothFR.traj*CPTOrth;
    CP = trajStruct(1).avgSmoothFR.traj*CPOrth;
    if CPT(end,3) < CPT(1,3)
        CPTOrth(:,3) = -1.*CPTOrth(:,3);
    end
    if CP(end,3) < CP(1,3)
        CPOrth(:,3) = -1.*CPOrth(:,3);
    end
    
    %Add Projections to trajStruct; do some stats 
    for i = 1:size(trajStruct,2)
        %Add average traces to trajStruct
        trajStruct(i).CIPCA.traj = (trajStruct(i).avgSmoothFR.traj-CIMargOffset')*CIPCA;
        trajStruct(i).targetPCA.traj = trajStruct(i).avgSmoothFR.traj*targetPCA;
        trajStruct(i).posturePCA.traj = trajStruct(i).avgSmoothFR.traj*posturePCA;
        trajStruct(i).postureLDA.traj = trajStruct(i).avgSmoothFR.traj*postureLDA;
        trajStruct(i).CPTOrth.traj = trajStruct(i).avgSmoothFR.traj*CPTOrth;
        trajStruct(i).CPOrth.traj = trajStruct(i).avgSmoothFR.traj*CPOrth;
        trajStruct(i).PTOrth.traj = trajStruct(i).avgSmoothFR.traj*PTOrth;
        trajStruct(i).CPTAllOrth.traj = trajStruct(i).avgSmoothFR.traj*CPTAllOrth;
    end
%%    
    %Plot CP and shadow 
    switch task
        case 'BCI'
            plotPostureList = [1:5];
            plotTargetList = [1];
            timePts = 1:minNumTimestamps;
            viewCoords = [-30,30];
    end
    figure
    hold on
     for posture = plotPostureList
       for target = plotTargetList
           if any([trajStruct.posture]==posture & [trajStruct.target]==target) 
               traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).CPOrth.traj(timePts,:); 
               plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
               plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
           end
       end
     end
     axis equal
     axis tight
     view([viewCoords])    
     ax = gca;
     zlim = ax.ZLim;
     zRange = zlim(2)-zlim(1);
     surfZCoord = zlim(1)-zRange/4; 
     for posture = plotPostureList
       for target = plotTargetList
           if any([trajStruct.posture]==posture & [trajStruct.target]==target) 
               traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).CPOrth.traj(timePts,:); 
               plot3(traj(:,1),traj(:,2),surfZCoord*ones(length(traj(:,1))),'Color',lightpcmap(posture,:),'LineWidth',2);
        
           end
       end
      end
     ax = gca; xlim = ax.XLim; ylim = ax.YLim;
     [X,Y] = meshgrid(xlim,ylim);
     Z = surfZCoord*ones(size(X,1),size(X,2));
     s = surf(X,Y,Z);
     s.FaceAlpha = 0.15; s.FaceColor = [0 0 1]; s.EdgeAlpha = 0;
    xticklabels({});yticklabels({});zticklabels({});
    xlabel(['Posture LDA'])
    zlabel(['Condition Invariant']) 
    grid on
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'cp.fig'));
        saveas(gcf,fullfile(saveDir,task,'cp.svg'));
    end
    xlabel([])
    zlabel([]) 
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'cpNoLabel.fig'));
        saveas(gcf,fullfile(saveDir,task,'cpNoLabel.svg'));
    end
    
%%     
    %Plot CPT with shadow 
    switch task
        case 'BCI'
            plotPostureList = [1:5];
            plotTargetList = [1,5];
            timePts = 1:minNumTimestamps;
            viewCoords = [-65,40];
    end
    figure
    hold on
     for posture = plotPostureList
       for target = plotTargetList
           if any([trajStruct.posture]==posture & [trajStruct.target]==target) 
               traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).CPTOrth.traj(timePts,:); 
               plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
               plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
           end
       end
     end
     axis equal
     axis tight
     view(viewCoords)
     ax = gca;
     xlim = ax.XLim; xRange = xlim(2)-xlim(1); xMid = mean(xlim); ax.XLim = [xMid-.75*xRange xMid+.75*xRange]; 
     ylim = ax.YLim; yRange = ylim(2)-ylim(1); yMid = mean(ylim); ax.YLim = [yMid-.75*yRange yMid+.75*yRange]; 
     zlim = ax.ZLim;
     zRange = zlim(2)-zlim(1);
     surfZCoord = zlim(1)-zRange/4;
      for posture = plotPostureList
       for target = plotTargetList
           if any([trajStruct.posture]==posture & [trajStruct.target]==target) 
               traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).CPTOrth.traj(timePts,:); 
               plot3(traj(:,1),traj(:,2),surfZCoord*ones(length(traj(:,1))),'Color',lightpcmap(posture,:),'LineWidth',2);
           end
       end
      end
     ax = gca; xlim = ax.XLim; ylim = ax.YLim; 
     [X,Y] = meshgrid(xlim,ylim);
     Z = surfZCoord*ones(size(X,1),size(X,2));
     s = surf(X,Y,Z);
     s.FaceAlpha = 0.15; s.FaceColor = [0 0 0]; s.EdgeAlpha = 0;
    xticklabels({}); yticklabels({}); zticklabels({});
    xlabel(['Posture LDA 1'])
    ylabel(['Target PC 1'])
    zlabel(['Condition Invariant']) 
    grid on
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'cpt.fig'));
        saveas(gcf,fullfile(saveDir,task,'cpt.svg'));
    end
    xlabel([]); ylabel([]); zlabel([]) 
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'cptNoLabel.fig'));
        saveas(gcf,fullfile(saveDir,task,'cptNoLabel.svg'));
    end
 %%   
    %Plot PT with shadow 
    switch task
        case 'BCI'
            plotPostureList = [1,3];
            plotTargetList = [1,3,5,7];
            timePts = 1:minNumTimestamps;
            viewCoords = [45,20];
    end
    figure
    hold on
     for posture = plotPostureList
       for target = plotTargetList
           if any([trajStruct.posture]==posture & [trajStruct.target]==target) 
               traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
               plot3(traj(:,2),traj(:,3),traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
               plot3(traj(1,2),traj(1,3),traj(1,1),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               plot3(traj(end,2),traj(end,3),traj(end,1),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
           end
       end
     end
     axis equal
     axis tight
     view(viewCoords)
     ax = gca;zlim = ax.ZLim;
     zRange = zlim(2)-zlim(1);
     surfZCoord = zlim(1)-zRange/3;
      for posture = plotPostureList
       for target = plotTargetList
           if any([trajStruct.posture]==posture & [trajStruct.target]==target) 
               traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
               plot3(traj(:,2),traj(:,3),surfZCoord*ones(length(traj(:,1))),'Color',lightpcmap(posture,:),'LineWidth',2);
           end
       end
      end
          ax = gca; xlim = ax.XLim; ylim = ax.YLim; 
    [X,Y] = meshgrid(xlim,ylim);
    Z = surfZCoord*ones(size(X,1),size(X,2));
    s = surf(X,Y,Z); s.FaceAlpha = 0.15; s.FaceColor = [1 0 0]; s.EdgeAlpha = 0;
    xticklabels({}); yticklabels({}); zticklabels({});
    xlabel(['Target PCA']); zlabel(['Posture LDA 1']) 
    grid on  
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'pt.fig'));
        saveas(gcf,fullfile(saveDir,task,'pt.svg'));
    end
    
    xlabel([]); zlabel([])  
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'ptNoLabel.fig'));
        saveas(gcf,fullfile(saveDir,task,'ptNoLabel.svg'));
    end
    %%
    
    %% Plot timecourses 
    switch task
        case 'BCI'
            plotPostureList = [1:5];
            plotTargetList = [3,7];
            timePts = 1:minNumTimestamps;
    end
    f = figure;
    [ha,pos] = tight_subplot(2,3,[.01 .03],[.1 .1],[.1 .01]) 
     f.Position = [50 50 900 300];
    hold on
     for posture = plotPostureList
       for target = plotTargetList
           if any([trajStruct.posture]==posture & [trajStruct.target]==target) 
               %traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).CPTAllOrth.traj(timePts,:); 
               time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(timePts); 
               axes(ha(1));
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).postureLDA.traj(timePts,:); 
                    plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                    hold on
               axes(ha(4));
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).postureLDA.traj(timePts,:); 
                    plot(time,traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                    hold on
               axes(ha(2));
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).targetPCA.traj(timePts,:); 
                    plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                    hold on
               axes(ha(5));
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).targetPCA.traj(timePts,:); 
                    plot(time,traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                    hold on
               axes(ha(3));
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).CIPCA.traj(timePts,:); 
                    plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                    hold on
               axes(ha(6));
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).CIPCA.traj(timePts,:); 
                    plot(time,traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                    hold on
           end
       end
     end
     set(ha(1:3),'XTickLabel','');
%      set(ha([2,3,5,6]),'YTickLabel','');
     axes(ha(1))
        title('Posture')
     axes(ha(2))
        title('Target')
     axes(ha(3))
        title('Condition Invariant')

%      subplot(4,1,1)
%         xticks({})
%         ylabel('Posture')
%      subplot(4,1,2)
%          xticks({})
%         ylabel('Target 1')
%      subplot(4,1,3)
%         xticks({})
%         ylabel('Target 2')
%      subplot(4,1,4)
%         xlabel('time (ms)')
%         ylabel('Condition-Invariant')
        
     %%
     %axis equal
     %axis tight
     
%     xticklabels({}); yticklabels({}); zticklabels({});
%     xlabel(['Target PCA']); zlabel(['Posture LDA 1']) 
%     grid on  
%     if saveFig
%         saveas(gcf,fullfile(saveDir,task,'pt.fig'));
%         saveas(gcf,fullfile(saveDir,task,'pt.svg'));
%     end
%     
%     xlabel([]); zlabel([])  
%     if saveFig
%         saveas(gcf,fullfile(saveDir,task,'ptNoLabel.fig'));
%         saveas(gcf,fullfile(saveDir,task,'ptNoLabel.svg'));
%     end

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
     
