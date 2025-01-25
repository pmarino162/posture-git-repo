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

%% Get posture and target lists
    postureList = unique([trajStruct.posture]);
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
    
%% Perform marginalizations and store them
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;

    
    %Do PCA on Marginalizations    
    allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
    Ca = cov(allTraj);
    [Ua,Sa] = eig(Ca);
    [Da,score,latent,tsquared,Expa,mu] = pca(allTraj);
    
    targetMargTraj = squeeze(targetMargTraj);
    targetMargTraj = reshape(targetMargTraj,[numTargets*minNumTimestamps,numChannels]);
    [Dt,score,latent,tsquared,Expt,mu] = pca(targetMargTraj); 
      
    postureMargTraj = squeeze(postureMargTraj);
    postureMargTraj = reshape(postureMargTraj,[numPostures*minNumTimestamps,numChannels]);
    [Dp,score,latent,tsquared,Expp,mu] = pca(postureMargTraj); 
    
    %Orthonormalize combinations of axes
    [PTOrth,~] = qr([Dp(:,1),Dt(:,1:2)]); PTOrth = PTOrth(:,1:3);
    
    
%% Top 3 PC projection
    figure
    hold on
    dims = [1,2,3];
    for posture = [1,5]
        for target = [1,3,5,7]
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            traj = traj(1:13,:);
            traj = traj*Da;
            plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
            plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
        end
    end

    axis equal
    view([-40 20])
    xticklabels({}); yticklabels({}); zticklabels({});
    grid on
    if saveFig
        saveas(gcf,fullfile(saveDir,'top3pc.fig'));
        saveas(gcf,fullfile(saveDir,'top3pc.svg'));
    end
    
%% Target Marginalization Embedded
    figure
    hold on
    for target = [1,3,5,7]
        ind = (target-1)*minNumTimestamps+1;
        traj = targetMargTraj(ind:ind+12,:);
        traj = traj*Dt;
        plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
        plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
        plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
    end
    axis equal
    axis tight
    view([-30 20])
    
    ax = gca;
    zlim = ax.ZLim;
    zRange = zlim(2)-zlim(1);
    surfZCoord = zlim(1)-zRange/4; 
    for target = [1,3,5,7]
        ind = (target-1)*minNumTimestamps+1;
        traj = targetMargTraj(ind:ind+12,:);
        traj = traj*Dt;
        plot3(traj(:,1),traj(:,2),surfZCoord*ones(length(traj(:,1))),'Color',lightpcmap(posture,:),'LineWidth',2);
    end
    ax = gca; xlim = ax.XLim; ylim = ax.YLim;
    [X,Y] = meshgrid(xlim,ylim);
    Z = surfZCoord*ones(size(X,1),size(X,2));
    s = surf(X,Y,Z);
    s.FaceAlpha = 0.15; s.FaceColor = [1 0 0]; s.EdgeAlpha = 0;
    xticklabels({}); yticklabels({}); zticklabels({});
    grid on
    if saveFig
        saveas(gcf,fullfile(saveDir,'targetMarg.fig'));
        saveas(gcf,fullfile(saveDir,'targetMarg.svg'));
    end
    
%% Posture Marginalization Embedded
    figure
    hold on
    for posture = [1:5]
        ind = (posture-1)*minNumTimestamps+1;
        traj = postureMargTraj(ind:ind+12,:);
        traj = traj*Dp;
        plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
        plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
        plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
    end
    axis equal
    axis tight
    view([-40 20])

    ax = gca;
    zlim = ax.ZLim;
    zRange = zlim(2)-zlim(1);
    surfZCoord = zlim(1)-zRange/4; 
    for posture = [1:5]
        ind = (posture-1)*minNumTimestamps+1;
        traj = postureMargTraj(ind:ind+12,:);
        traj = traj*Dp;
        plot3(traj(:,1),traj(:,2),surfZCoord*ones(length(traj(:,1))),'Color',lightpcmap(posture,:),'LineWidth',2);
    end
    ax = gca; xlim = ax.XLim; ylim = ax.YLim;
    [X,Y] = meshgrid(xlim,ylim);
    Z = surfZCoord*ones(size(X,1),size(X,2));
    s = surf(X,Y,Z);
    s.FaceAlpha = 0.15; s.FaceColor = [0 0 1]; s.EdgeAlpha = 0;
    xticklabels({}); yticklabels({}); zticklabels({});
    grid on
    if saveFig
        saveas(gcf,fullfile(saveDir,'postureMarg.fig'));
        saveas(gcf,fullfile(saveDir,'postureMarg.svg'));
   end
    
%% Target Marginalizaiton in Target Subspace
    figure
    hold on
    for target = [1,3,5,7]
        ind = (target-1)*minNumTimestamps+1;
        traj = targetMargTraj(ind:ind+12,:);
        traj = traj*Dt;
        plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
        plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
        plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
    end
    axis equal
    view([0 90])
    ax = gca; txlim = ax.XLim; tylim = ax.YLim; tzlim = ax.ZLim;
    [Xt,Yt] = meshgrid(txlim,tylim);
    surfZCoord = tzlim(1); 
    Zt = surfZCoord*ones(size(Xt,1),size(Xt,2));
    s = surf(Xt,Yt,Zt);     
    s.FaceAlpha = 0.15; s.FaceColor = [1 0 0]; s.EdgeAlpha = 0;
    xticklabels({}); yticklabels({}); zticklabels({});
    grid on
    if saveFig
        saveas(gcf,fullfile(saveDir,'targTargSubspace.fig'));
        saveas(gcf,fullfile(saveDir,'targTargSubspace.svg'));
    end
   
%% Posture Marginalizaiton in Posture Subspace
    figure
    hold on
    for posture = [1:5]
        ind = (posture-1)*minNumTimestamps+1;
        traj = postureMargTraj(ind:ind+12,:);
        traj = traj*Dp;
        plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
        plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
        plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
    end
    axis equal
    view([0 90])
    ax = gca; pxlim = ax.XLim; pylim = ax.YLim; pzlim = ax.ZLim;
    [Xp,Yp] = meshgrid(pxlim,pylim);
    surfZCoord = pzlim(1); 
    Zp = surfZCoord*ones(size(Xp,1),size(Xp,2));
    s = surf(Xt,Yt,Zt);
    s.FaceAlpha = 0.15; s.FaceColor = [0 0 1]; s.EdgeAlpha = 0;
    xticklabels({}); yticklabels({}); zticklabels({});
    grid on    
    if saveFig
        saveas(gcf,fullfile(saveDir,'postPostSubspace.fig'));
        saveas(gcf,fullfile(saveDir,'postPostSubspace.svg'));
    end
    
%% Target Marginalizaiton in Posture Subspace
    figure
    hold on
    for target = [1,3,5,7]
        ind = (target-1)*minNumTimestamps+1;
        traj = targetMargTraj(ind:ind+12,:);
        traj = traj*Dp;
        plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
        plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
        plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
    end
    
    axis equal
    view([0 90])
    s = surf(Xt,Yt,Zt);
    s.FaceAlpha = 0.15; s.FaceColor = [0 0 1]; s.EdgeAlpha = 0;
    xticklabels({}); yticklabels({}); zticklabels({});
    grid on
    if saveFig
        saveas(gcf,fullfile(saveDir,'targPostSubspace.fig'));
        saveas(gcf,fullfile(saveDir,'targPostSubspace.svg'));
    end
    
%% Posture Marginalizaiton in Target Subspace
    figure
    hold on
    for posture = [1:5]
        ind = (posture-1)*minNumTimestamps+1;
        traj = postureMargTraj(ind:ind+12,:);
        traj = traj*Dt;
        plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
        plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
        plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
    end
    axis equal
    view([0 90])
    s = surf(Xt,Yt,Zt);
    s.FaceAlpha = 0.15; s.FaceColor = [1 0 0]; s.EdgeAlpha = 0;
    xticklabels({}); yticklabels({}); zticklabels({});
    grid on
    if saveFig
        saveas(gcf,fullfile(saveDir,'postTargSubspace.fig'));
        saveas(gcf,fullfile(saveDir,'postTargSubspace.svg'));
    end
