clear; clc; clf; close all;
 
%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\BCI Figures';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;
    pcmap = parula(5);
    lightpcmap = rgb2hsv(pcmap);
    lightpcmap(:,2)=.3;
    lightpcmap =hsv2rgb(lightpcmap);
    
%% Load Data, subselect, get traj struct
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;
    
    [Data,I30GPFAParams,I30DecoderParams,...
     N00GPFAParams,N00DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200314_20211210;

     condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    
    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
        trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',0},{'state','Step 2','first',0}};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(2).inclStates = {{'state','Reach','first',0},{'state','Hold','first',0}};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
    
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
 
%% Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);
    
%% Peform LDA by posture on BCI and iso force
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = postureList
       obsStruct(structInd).label = posture;
       allObs = [];
       tempTrajStruct = trajStruct([trajStruct.posture]==posture);
       tempTrajStruct = tempTrajStruct([tempTrajStruct.task]==1 | [tempTrajStruct.task]==3);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
               traj = tempTrajStruct(i).allSmoothFR(j).traj;
%                traj = traj(1:minNumTimestamps,:);
               allObs = vertcat(allObs,traj);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    
    [postureLDA] = doLDA(obsStruct);   
    
%%    one task at a time
for task = 1:3    

    tempTrajStruct = trajStruct([trajStruct.task]==task);
    targetList = unique([tempTrajStruct.target]); 
    numTargets = size(targetList,2);
    % Get timestamps dist and min number timestamps
    numTimestamps = [];
    numCondTraj = [];
    for i = 1:size(tempTrajStruct,2)
       numTraj = size(tempTrajStruct(i).allSmoothFR,2);
       numCondTraj = [numCondTraj,numTraj];
       for j = 1:numTraj
          numTimestamps = [numTimestamps,size(tempTrajStruct(i).allSmoothFR(j).timestamps,2)]; 
       end
    end
    figure
    histogram(numTimestamps)
    xlabel('Number of 25ms bins')
    ylabel('Number of trials')
    
    figure
    histogram(numCondTraj)
    xlabel('Number of trials')
    ylabel('Number of conditions')
    
    % Get minimum number of trials and timestamps
    [minNumTimestamps,i] = min(numTimestamps);
    [minNumCondTraj,i] = min(numCondTraj);
    [maxNumCondTrials,i] = max(numCondTraj);
    

    

% Form X, containing trial-averaged data for each condition
    X = NaN(minNumTimestamps,numTargets,numPostures,maxNumCondTrials,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            numTrials = size(tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).allSmoothFR,2);
            for trial = 1:numTrials
                traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).allSmoothFR(trial).traj;
                X(:,targetInd,postureInd,trial,:) = traj(1:minNumTimestamps,:); 
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    X = squeeze(mean(X,4,'omitnan'));
    
% Perform marginalizations and store them
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
    
% Different way to do marginalizations
    targetTrajNoCP = X - CIMargTraj - CIMargOffset - postureMargTraj - postureMargOffset;
    
% Do PCA on Marginalizations    
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
    

    
% Orthonormalize combinations of axes
    CPTOrth = qr([postureLDA(:,1),targetNoCPPCA(:,1),CIPCA(:,1)]);
%     CPTOrth = qr([postureLDA(:,1),targetNoCPPCA(:,1:2)]);
    CPTOrth = CPTOrth(:,1:3);
    
    CPOrth = qr([postureLDA(:,1:2),CIPCA(:,1)]);
    CPOrth = CPOrth(:,1:3);
    PTOrth = qr([postureLDA(:,1),targetPCA(:,1:2)]);
    PTOrth = PTOrth(:,1:3);
    
% Add Projections to trajStruct; do some stats 
    for i = 1:size(tempTrajStruct,2)
        %Add average traces to trajStruct
        tempTrajStruct(i).CIPCA.traj = (tempTrajStruct(i).avgSmoothFR.traj-CIMargOffset')*CIPCA;
        tempTrajStruct(i).targetPCA.traj = tempTrajStruct(i).avgSmoothFR.traj*targetPCA;
        tempTrajStruct(i).posturePCA.traj = tempTrajStruct(i).avgSmoothFR.traj*posturePCA;
        tempTrajStruct(i).CPTOrth.traj = tempTrajStruct(i).avgSmoothFR.traj*CPTOrth;
        tempTrajStruct(i).CPOrth.traj = tempTrajStruct(i).avgSmoothFR.traj*CPOrth;
        tempTrajStruct(i).PTOrth.traj = tempTrajStruct(i).avgSmoothFR.traj*PTOrth;
    end
    

% Plot CP and shadow 
    switch task
        case 1
            plotPostureList = [1:5];
            plotTargetList = [1];
            timePts = 1:20;
        case 2
            plotPostureList = [1:5];
            plotTargetList = [1];
            timePts = 1:10;
        case 3
            plotPostureList = [1:5];
            plotTargetList = [1];
            timePts = 1:10;
    end

    figure
   hold on
     for posture = plotPostureList
       for target = plotTargetList
           if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target) 
               traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).CPOrth.traj(timePts,:); 
               plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
               plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
           end
       end
     end
     
     ax = gca;
     
     view([-45,20])
     
     ax = gca;
     xlim = ax.XLim;
     ylim = ax.YLim;
     zlim = ax.ZLim;
     
     surfZCoord = zlim(1)-100;
   
     [X,Y] = meshgrid(xlim,ylim);
     Z = surfZCoord*ones(size(X,1),size(X,2));
     s = surf(X,Y,Z);
     s.FaceAlpha = 0.15;
     s.FaceColor = [0 0 1];
     s.EdgeAlpha = 0;
     
      for posture = plotPostureList
       for target = plotTargetList
           if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target) 
               traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).CPOrth.traj(timePts,:); 
               plot3(traj(:,1),traj(:,2),surfZCoord*ones(length(traj(:,1))),'Color',lightpcmap(posture,:),'LineWidth',2);
        
           end
       end
      end
     
         xticklabels({});
    yticklabels({});
    zticklabels({});
    xlabel(['Posture LDA'])
    zlabel(['Condition Invariant']) 
    grid on
%     axis equal


% Plot CPT with shadow 
    plotPostureList = [1:5];
    plotTargetList = [1,5];
    timePts = 1:12;

    figure
   hold on
     for posture = plotPostureList
       for target = plotTargetList
           if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target) 
               traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).CPTOrth.traj(timePts,:); 
               plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
               plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
           end
       end
     end
     
     view([-45,40])
     
     ax = gca;
     xlim = ax.XLim;
     ylim = ax.YLim;
     zlim = ax.ZLim;
     
     surfZCoord = zlim(1)-200;
   
     [X,Y] = meshgrid(xlim,ylim);
     Z = surfZCoord*ones(size(X,1),size(X,2));
     s = surf(X,Y,Z);
     s.FaceAlpha = 0.15;
     s.FaceColor = [0 0 0];
     s.EdgeAlpha = 0;
    
      for posture = plotPostureList
       for target = plotTargetList
           if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target) 
               traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).CPTOrth.traj(timePts,:); 
               plot3(traj(:,1),traj(:,2),surfZCoord*ones(length(traj(:,1))),'Color',lightpcmap(posture,:),'LineWidth',2);
        
           end
       end
      end

    xticklabels({});
    yticklabels({});
    zticklabels({});
    xlabel(['Posture LDA 1'])
    ylabel(['Target PC 1'])
    zlabel(['Condition Invariant']) 
    grid on
    

    
% Plot PT with shadow 
    plotPostureList = [1,5];
    plotTargetList = [1,3,5,7];
    timePts = 1:11;

    figure
    hold on
     for posture = plotPostureList
       for target = plotTargetList
           if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target) 
               traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth.traj(timePts,:); 
               plot3(traj(:,2),traj(:,3),traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
               plot3(traj(1,2),traj(1,3),traj(1,1),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               plot3(traj(end,2),traj(end,3),traj(end,1),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
           end
       end
     end
     
     view([45,20])
     
     ax = gca;
     xlim = ax.XLim;
     ylim = ax.YLim;
     zlim = ax.ZLim;
     
     surfZCoord = zlim(1)-100;
   
     [X,Y] = meshgrid(xlim,ylim);
     Z = surfZCoord*ones(size(X,1),size(X,2));
     s = surf(X,Y,Z);
     s.FaceAlpha = 0.15;
     s.FaceColor = [1 0 0];
     s.EdgeAlpha = 0;
     
      for posture = plotPostureList
       for target = plotTargetList
           if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target) 
               traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth.traj(timePts,:); 
               plot3(traj(:,2),traj(:,3),surfZCoord*ones(length(traj(:,1))),'Color',lightpcmap(posture,:),'LineWidth',2);
           end
       end
     end
     
    xticklabels({});
    yticklabels({});
    zticklabels({});
    xlabel(['Target PCA'])
    zlabel(['Posture LDA 1']) 
    grid on  

end

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
     
