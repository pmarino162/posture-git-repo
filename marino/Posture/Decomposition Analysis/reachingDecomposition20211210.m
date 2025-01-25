clear; clc; clf; close all;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;
    
%% Load Data 
    [Data] = loadEarlData20210706_20211210;
    
%% Get trajStruct 
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'GridReaching'};
    binWidth = 25;

    %Reach
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-300},{'state','Target Hold','first',400}};
    reachTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    %Hold
    trialInclStates(1).inclStates = {{'state','Target Hold','first',150},{'state','Success with Reward','first',0}};
    holdTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%% Get posture and target lists
    trajStruct = reachTrajStruct;
    postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
   
%% Perform posture LDA on holding neural activity numDims
    trajStruct = holdTrajStruct([holdTrajStruct.target]==3);
    trajStruct = trajStruct(ismember([trajStruct.posture],[1:4]));
    %Get number of points for each condition
    obsStruct = struct('target',[],'posture',[],'numObs',[],'allObs',[]);
    structInd = 1;
    for i = 1:size(trajStruct,2)
       obsStruct(structInd).target = trajStruct(i).target;
       obsStruct(structInd).posture = trajStruct(i).posture;
       obsStruct(structInd).allObs = vertcat(trajStruct(i).allSmoothFR.traj);
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    
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
        labels(k:k+minNumObs-1,1) = ones(minNumObs,1).*obsStruct(i).posture;
        k = k+minNumObs;
    end

    postureLDA = fisherLDA(obs, labels);
    [postureLDA,~] = qr(postureLDA);
    postureLDA = postureLDA(:,1:3);
    
% %% Try posture LDA on neural activity at movement onset
%     trajStruct = reachTrajStruct(ismember([reachTrajStruct.posture],[1:4]));
%     trajStruct = trajStruct(ismember([trajStruct.target],[3]));
%     %Get number of points for each condition
%     obsStruct = struct('posture',[],'numObs',[],'allObs',[]);
%     structInd = 1;
%     for posture = [1:4]
%        tempTrajStruct = trajStruct([trajStruct.posture]==posture);
%        obsStruct(structInd).posture = posture;
%        allObsInd = 1;
%        for i = 1:size(tempTrajStruct,2)
%            numTraj = size(tempTrajStruct(i).allSmoothFR,2);
%            for j = 1:numTraj
%                zeroInd = find(tempTrajStruct(i).allSmoothFR(j).timestamps == 0);
%                obsStruct(structInd).allObs(allObsInd:allObsInd+4,:) = tempTrajStruct(i).allSmoothFR(j).traj(zeroInd-2:zeroInd+2,:);
%                allObsInd = allObsInd + 5;
%            end
%        end
%        obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
%        structInd = structInd + 1;
%     end
%     
%     %Preallocate
%     minNumObs = min([obsStruct.numObs]);
%     numClasses = size(obsStruct,2);
%     numDims = size(obsStruct(1).allObs,2);
%     obs = NaN(minNumObs*numClasses,numDims); labels =  NaN(minNumObs*numClasses,1); 
% 
%     %Fill
%     k = 1;
%     for i = 1:size(obsStruct,2)
%         totalClassObs = size(obsStruct(i).allObs,1);
%         obs(k:k+minNumObs-1,:) = obsStruct(i).allObs(randsample(totalClassObs,minNumObs),:);
%         labels(k:k+minNumObs-1,1) = ones(minNumObs,1).*obsStruct(i).posture;
%         k = k+minNumObs;
%     end
% 
%     postureLDA = fisherLDA(obs, labels);
%     [postureLDA,~] = qr(postureLDA);
%     postureLDA = postureLDA(:,1:3);
    
%% Peform target LDA on neural activity at movement onset
%     trajStruct = reachTrajStruct([reachTrajStruct.posture]==2);%reachTrajStruct;%([reachTrajStruct.posture]==5);
    trajStruct = reachTrajStruct;
    %Get number of points for each condition
    reachObsStruct = struct('target',[],'numObs',[],'allObs',[]);
    structInd = 1;
    for target = targetList
       tempTrajStruct = trajStruct([trajStruct.target]==target);
       reachObsStruct(structInd).target = target;
       allObsInd = 1;
       for i = 1:size(tempTrajStruct,2)
           numTraj = size(tempTrajStruct(i).allSmoothFR,2);
           for j = 1:numTraj
               zeroInd = find(tempTrajStruct(i).allSmoothFR(j).timestamps == 0);
               reachObsStruct(structInd).allObs(allObsInd:allObsInd+4,:) = tempTrajStruct(i).allSmoothFR(j).traj(zeroInd-2:zeroInd+2,:);
               allObsInd = allObsInd + 5;
           end
       end
       reachObsStruct(structInd).numObs = size(reachObsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    
    %Preallocate
    minNumObs = min([reachObsStruct.numObs]);
    numClasses = size(reachObsStruct,2);
    numDims = size(reachObsStruct(1).allObs,2);
    reachObs = NaN(minNumObs*numClasses,numDims); reachLabels =  NaN(minNumObs*numClasses,1); 

    %Fill
    k = 1;
    for i = 1:size(reachObsStruct,2)
        totalClassObs = size(reachObsStruct(i).allObs,1);
        reachObs(k:k+minNumObs-1,:) = reachObsStruct(i).allObs(randsample(totalClassObs,minNumObs),:);
        reachLabels(k:k+minNumObs-1,1) = ones(minNumObs,1).*reachObsStruct(i).target;
        k = k+minNumObs;
    end

    targetLDA = fisherLDA(reachObs, reachLabels);
    [targetLDA,~] = qr(targetLDA);
    targetLDA = targetLDA(:,1:7);

%% Do PCA on Condition-averaged reaching data
    avgSmoothFR = [reachTrajStruct.avgSmoothFR];
    allAvgs = vertcat(avgSmoothFR.traj);
    [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);
    
%% Do PCA in null space of posture space 
    postureAxisNullBasis = null(postureLDA(:,1:2)');
    nullProj = allAvgs*postureAxisNullBasis;
    [postureAxisNullBasisPC,score,latent,tsquared,explained,mu] = pca(nullProj);
    postureAxisNullBasis = postureAxisNullBasis*postureAxisNullBasisPC;
    
%% Get Post Targ Orth
    [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:3);
    
%% Get angles
    acosd(dot(postureLDA(:,1),targetLDA(:,1)))
    acosd(dot(postureLDA(:,1),targetLDA(:,2)))
    
%% Add Projections to trajStruct and modelStruct
    for i = 1:size(reachTrajStruct,2)
       reachTrajStruct(i).PCA = reachTrajStruct(i).avgSmoothFR.traj*coeff;
       reachTrajStruct(i).postureLDA = reachTrajStruct(i).avgSmoothFR.traj*postureLDA;
       reachTrajStruct(i).targetLDA = reachTrajStruct(i).avgSmoothFR.traj*targetLDA;
       reachTrajStruct(i).PCApNull = reachTrajStruct(i).avgSmoothFR.traj*postureAxisNullBasis;
       reachTrajStruct(i).PTOrth = reachTrajStruct(i).avgSmoothFR.traj*postTargOrth;
%        holdTrajStruct(i).postureLDA = holdTrajStruct(i).avgSmoothFR.traj*postureLDA;
%        holdTrajStruct(i).targetLDA = holdTrajStruct(i).avgBinFR.traj*targetLDA;
    end    
    
%% Plot Posture Hold points and reach traj
     figure
    for i = 1:size(obs,1)
        pt = obs(i,:)*postureLDA;
        label = labels(i);
        plot3(pt(:,1),pt(:,2),pt(:,3),'.','MarkerSize',10,'Color',pcmap(label,:));
        hold on
    end
    xlabel('Posture LDA 1')
    ylabel('Posture LDA 2')
    zlabel('Posture LDA 3')
    axis equal
    
    figure
    for posture = 2
       for target = targetList
          traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgSmoothFR.traj*postureLDA; 
          plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(target,:));
          hold on
          plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',15,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
          plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',5,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
       end
    end
    xlabel('Posture LDA 1')
    ylabel('Posture LDA 2')
    zlabel('Posture LDA 3')
    
%% Plot Target Training points and reach traj
    figure
    for i = 1:size(reachObs,1)
        pt = reachObs(i,:)*targetLDA;
        label = reachLabels(i);
        plot3(pt(:,1),pt(:,2),pt(:,3),'.','MarkerSize',10,'Color',tcmap(label,:));
        hold on
    end
    xlabel('Target LDA 1')
    ylabel('Target LDA 2')
    zlabel('Target LDA 3')
    
    figure
    for posture = 2
       for target = targetList
          traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgSmoothFR.traj*targetLDA; 
          plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(target,:));
          hold on
          plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',15,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
          plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',5,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
       end
    end
    xlabel('Target LDA 1')
    ylabel('Target LDA 2')
    zlabel('Target LDA 3')
    
%% Plot PCA Trajectories
    for posture = 1:4
       for target = 1:8
           traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).PCA; 
           plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(target,:));
           hold on
       end
    end
    xlabel('PC 1')
    ylabel('PC 2')
    zlabel('PC 3')
 
%% Plot PCA Timecourses
    for target = 1:8
        figure
        hold on
        for posture = 1:7
            if ~isempty(reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target))
               time = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgSmoothFR.timestamps; 
               postureTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).postureLDA; 
               pNullTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).PCApNull; 
               traj = [postureTraj(:,1),pNullTraj];
               for i = 1:10
                   subplot(2,5,i)
                   plot(time,traj(:,i),'Color',pcmap(posture,:),'LineWidth',1.5);
                   hold on
               end
               sgtitle(['Target ',num2str(target)])
            end
        end
        
            %Get axis ranges
    for plotInd = 1:10
       ax = subplot(2,5,plotInd);
       range(plotInd) = ax.YLim(1,2)-ax.YLim(1,1);
       center(plotInd) = (ax.YLim(1,2)+ax.YLim(1,1))/2;
    end
    maxRange = max(range);
    for plotInd = 1:10
       ax = subplot(2,5,plotInd);
       ax.YLim = [center(plotInd)-maxRange/2 center(plotInd)+maxRange/2];
       xlabel('time (ms)')
       if plotInd == 1
           ylabel('Posture LDA 1')
       else
        ylabel(['PC ',num2str(plotInd-1)])
       end
    end
        
        
        
    end

     

    
    
%% Posture 1 and top 2 PCs
    for target = [1:8]
        figure
        hold on
        for posture = 1:4
           postureTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).postureLDA; 
           PCTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).PCApNull; 
           plot3(postureTraj(:,1),PCTraj(:,1),PCTraj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
           plot3(postureTraj(1,1),PCTraj(1,1),PCTraj(1,2),'.','MarkerSize',15,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
           plot3(postureTraj(end,1),PCTraj(end,1),PCTraj(end,2),'>','MarkerSize',5,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
        end
        xlabel('Posture LDA 1')
        ylabel('PC 1')
        zlabel('PC 2')
        title(['Target ',num2str(target)])
    end


%% PxT Projection 
 
    figure
    hold on
     for posture = 2
       for target = [1:8]
           traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).PTOrth; 
           plot3(traj(:,1),traj(:,2),traj(:,3),'Color','k','LineWidth',2);
           plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
           plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',10,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
       end
    end
    xlabel('Posture LDA 1')
    ylabel('Target LDA 1')
    zlabel('Target LDA 2')
    
%% PxT Time course 
    figure
    hold on
 
     for posture = [1,4]
            targInd = 1;
       for target = [2,6]
           time = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgSmoothFR.timestamps; 
           traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).PTOrth; 
           for i = 1:3
               subplot(1,3,i)
               if targInd ==1 
                plot(time,traj(:,i),'Color',tcmap(posture,:),'LineWidth',1.5);
               else
                   plot(time,traj(:,i),'--','Color',tcmap(posture,:),'LineWidth',1.5);
               end
               hold on
           end
           targInd = targInd + 1;
       end
     end

     for i = 1:3
         ax = subplot(1,3,i);
         xlabel('time (ms)')
         if i==1
             ylabel('Posture LDA')
         elseif i==2
             ylabel('Target LDA 1')
         elseif i==3
             ylabel('Target LDA 2')
         end
     end