clear; clc; clf; close all;
    
%% Set up figure saving
    saveFig = false; 
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220106 - planning decomposition';
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = customRainbow;
    
%% Load Data 
    [Data] = loadEarlData20210706_20211210;
    
%% Get trajStruct 
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'GridReaching'};
    binWidth = 25;

    %Delay
    kinData = [Data.kinData];
    delayLength = [kinData.delayLength];
    delayData = Data(delayLength>=500);
    trialInclStates(1).inclStates = {{'state','Delay','first',0},{'state','Target Acquire','first',0}};
    delayTrajStruct = getTrajStruct20211210(delayData,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    %Reach
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-300},{'state','Target Hold','first',400}};
    reachTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    %Hold
    trialInclStates(1).inclStates = {{'state','Target Hold','first',150},{'state','Success with Reward','first',0}};
    holdTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%% Do PCA on Condition-averaged delay and reaching data; add to trajStructs
    %Delay
    avgSmoothFR = [delayTrajStruct.avgSmoothFR];
    allDelayAvgs = vertcat(avgSmoothFR.traj);
    [delayCoeff,~,~,~,delayExplained,delayMu] = pca(allDelayAvgs);
    for i = 1:size(delayTrajStruct,2)
        delayTrajStruct(i).delayPCA = delayTrajStruct(i).avgSmoothFR.traj*delayCoeff;
    end

    %Reach
    avgSmoothFR = [reachTrajStruct.avgSmoothFR];
    allReachAvgs = vertcat(avgSmoothFR.traj);
    [reachCoeff,~,~,~,reachExplained,reachMu] = pca(allReachAvgs);
    for i = 1:size(reachTrajStruct,2)
       reachTrajStruct(i).reachPCA = reachTrajStruct(i).avgSmoothFR.traj*reachCoeff;
    end   
    
%% Plot PCA Trajectories
    %All targets for each posture
    plotTarget = 1:8;
    for plotPosture = 1:7
        %3D PCA Plot
        figure
        for posture = plotPosture
           for target = plotTarget
               if any([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target)
                   traj = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).delayPCA; 
                   plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(target,:));
                   hold on
                   plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',15,'Color',tcmap(target,:));
                   plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',5,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',tcmap(target,:));
               end
           end
        end
        xlabel('Delay PC 1')
        ylabel('Delay PC 2')
        zlabel('Delay PC 3')   
        grid on
        title(['Posture ',num2str(plotPosture)])
        axis equal
        if saveFig
            saveas(gcf,[saveDir,'\Post',num2str(plotPosture),'PCA3.fig']);
        end
        
        f = figure; f.Position = [0 0 1500 500];
        for posture = plotPosture
            for target = plotTarget
                if any([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target)
                    traj = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).delayPCA; 
                    time = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).avgSmoothFR.timestamps; 
                    for i = 1:10
                       subplot(2,5,i)
                       plot(time,traj(:,i),'Color',tcmap(target,:),'LineWidth',1);
                       hold on
                    end
                end
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
           ylabel(['PC ',num2str(plotInd)])
        end
        sgtitle(['Posture ',num2str(plotPosture)])
        if saveFig
            saveas(gcf,[saveDir,'\Post',num2str(plotPosture),'PCA.fig']);
            saveas(gcf,[saveDir,'\Post',num2str(plotPosture),'PCA.jpg']);
        end
    end
    
    %All Targets for each Posture
    plotPosture = 1:7;
    for plotTarget = 1:8
        %3D PCA Plot
        figure
        for posture = plotPosture
           for target = plotTarget
               if any([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target)
                   traj = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).delayPCA; 
                   plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(posture,:));
                   hold on
                   plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',15,'Color',tcmap(posture,:));
                   plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',5,'MarkerFaceColor',tcmap(posture,:),'MarkerEdgeColor',tcmap(posture,:));
               end
           end
        end
        xlabel('Delay PC 1')
        ylabel('Delay PC 2')
        zlabel('Delay PC 3')   
        grid on
        title(['Target ',num2str(plotTarget)])
        axis equal
        if saveFig
            saveas(gcf,[saveDir,'\Targ',num2str(plotTarget),'PCA3.fig']);
        end
        
        f = figure; f.Position = [0 0 1500 500];
        for posture = plotPosture
            for target = plotTarget
                if any([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target)
                    traj = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).delayPCA; 
                    time = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).avgSmoothFR.timestamps; 
                    for i = 1:10
                       subplot(2,5,i)
                       plot(time,traj(:,i),'Color',tcmap(posture,:),'LineWidth',1);
                       hold on
                    end
                end
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
           ylabel(['PC ',num2str(plotInd)])
        end
        sgtitle(['Target ',num2str(plotTarget)])
        if saveFig
            saveas(gcf,[saveDir,'\Targ',num2str(plotTarget),'PCA.fig']);
            saveas(gcf,[saveDir,'\Targ',num2str(plotTarget),'PCA.jpg']);
        end
    end
    
%     %Delay
%     figure
%     for posture = [1,4]
%        for target = [2,4,6,8]
%            if any([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target)
%                traj = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).delayPCA; 
%                plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(posture,:));
%                hold on
%                plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',15,'Color',tcmap(target,:));
%                plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',5,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',tcmap(target,:));
%            end
%        end
%     end
%     xlabel('Delay PC 1')
%     ylabel('Delay PC 2')
%     zlabel('Delay PC 3')   
%     
%     %Reach
%     figure
%     for posture = [1,4]
%        for target = [2,4,6,8]
%            if any([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target)
%                traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).reachPCA; 
%                plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(target,:));
%                hold on
%                plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',15,'Color',tcmap(target,:));
%                plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',5,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',tcmap(target,:));
%            end
%        end
%     end
%     xlabel('Reach PC 1')
%     ylabel('Reach PC 2')
%     zlabel('Reach PC 3')    
%     

%% "Subtract IC" plot
    plotPosture = [1,7];
    plotTarget = 1:8;
    figure
    for posture = plotPosture
       for target = plotTarget
           if any([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target)
               traj = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).delayPCA; 
               traj = traj - traj(1,:);
               plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(target,:));
               hold on
               plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',15,'Color',tcmap(target,:));
               plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',5,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',tcmap(target,:));
           end
       end
    end
    xlabel('Delay PC 1')
    ylabel('Delay PC 2')
    zlabel('Delay PC 3')   
    grid on
    title(['Subtract ICs'])
    axis equal
    if saveFig
        saveas(gcf,[saveDir,'\SubtractICs_PCA3.fig']);
    end
        
%% Get posture and target lists
    trajStruct = reachTrajStruct;
    postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
   
%% Peform LDA by posture on all data 
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = 1:4
       obsStruct(structInd).label = posture;
       allObs = [];
       tempTrajStruct = delayTrajStruct([delayTrajStruct.posture]==posture);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
               traj = tempTrajStruct(i).allSmoothFR(j).traj;
               traj = traj(1:20,:);
               allObs = vertcat(allObs,traj);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    
    [postureLDA] = doLDA(obsStruct);

  
%% Peform LDA by target on all data - around 250ms after cue onset
    trajStruct = delayTrajStruct(ismember([delayTrajStruct.posture],[1:4]));
    %Get number of points for each condition
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;
    for target = 1:8
       obsStruct(structInd).label = target;
       allObs = [];
       tempTrajStruct = delayTrajStruct([delayTrajStruct.target]==target);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
               traj = tempTrajStruct(i).allSmoothFR(j).traj;
               [~,ind275] = min(abs(timestamps-275));
               traj = traj(ind275-2:ind275+2,:);
               allObs = vertcat(allObs,traj);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    
    [targetLDA] = doLDA(obsStruct);
    
    
%% Do PCA in null space of posture space 
    postureAxisNullBasis = null(postureLDA(:,1:2)');
    nullProj = allDelayAvgs*postureAxisNullBasis;
    [postureAxisNullBasisPC,score,latent,tsquared,explained,mu] = pca(nullProj);
    postureAxisNullBasis = postureAxisNullBasis*postureAxisNullBasisPC;
    
%% Get Post Targ Orth
    [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:3);
    
%% Get angles
    acosd(dot(postureLDA(:,1),targetLDA(:,1)))
    acosd(dot(postureLDA(:,1),targetLDA(:,2)))
    
%% Add Projections to trajStruct 
    %Delay
    for i = 1:size(delayTrajStruct,2)
        delayTrajStruct(i).postureLDA = delayTrajStruct(i).avgSmoothFR.traj*postureLDA;
        delayTrajStruct(i).targetLDA = delayTrajStruct(i).avgSmoothFR.traj*targetLDA;
        delayTrajStruct(i).PCApNull = delayTrajStruct(i).avgSmoothFR.traj*postureAxisNullBasis;
        delayTrajStruct(i).PTOrth = delayTrajStruct(i).avgSmoothFR.traj*postTargOrth;
    end

    %Hold
    for i = 1:size(holdTrajStruct,2)
%        holdTrajStruct(i).postureLDA = holdTrajStruct(i).avgSmoothFR.traj*postureLDA;
%        holdTrajStruct(i).targetLDA = holdTrajStruct(i).avgBinFR.traj*targetLDA;
    end    
    
%% Do plots involving LDA Projections 
    % Posture LDA and top 2 PCs
    figure
    hold on
     for posture = 1:4
       for target = [1,5]
           postureLDA = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).postureLDA; 
           PCApNull = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).PCApNull;
           traj = [postureLDA(:,1),PCApNull(:,1:2)];
           plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(posture,:),'LineWidth',1);
           plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
           plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',10,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
       end
    end
    xlabel('Posture LDA 1')
    ylabel('Delay PC 1')
    zlabel('Delay PC 2') 
     axis equal

   % PxT Projection 
    figure
    hold on
     for posture = 1:4
       for target = 1:8
           traj = delayTrajStruct([delayTrajStruct.posture]==posture & [delayTrajStruct.target]==target).PTOrth; 
           plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(posture,:),'LineWidth',1.5);
           plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',10,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
           plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',5,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
       end
    end
    xlabel('Posture LDA 1')
    ylabel('Target LDA 1')
    zlabel('Target LDA 2')
    axis equal
    
    

%% Plot Posture Hold points and reach traj
     figure
    for i = 1:size(obs,1)
        pt = obs(i,:)*postureLDA;
        label = labels(i);
        plot3(pt(:,1),pt(:,2),pt(:,3),'.','MarkerSize',10,'Color',tcmap(label,:));
        hold on
    end
    xlabel('Posture LDA 1')
    ylabel('Posture LDA 2')
    zlabel('Posture LDA 3')
    
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
                   plot(time,traj(:,i),'Color',tcmap(posture,:),'LineWidth',1);
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
    for target = [1,3,5]
        figure
        hold on
        for posture = 1:7
           postureTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).postureLDA; 
           PCTraj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).PCApNull; 
           plot3(postureTraj(:,1),PCTraj(:,1),PCTraj(:,2),'Color',tcmap(posture,:),'LineWidth',2);
           plot3(postureTraj(1,1),PCTraj(1,1),PCTraj(1,2),'.','MarkerSize',15,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
           plot3(postureTraj(end,1),PCTraj(end,1),PCTraj(end,2),'>','MarkerSize',5,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
        end
        xlabel('Posture LDA 1')
        ylabel('PC 1')
        zlabel('PC 2')
        title(['Target ',num2str(target)])
    end



    
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