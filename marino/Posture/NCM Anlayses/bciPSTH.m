clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\BCI_PSTH';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Run loop for each dataset    
    for datasetList = {'E20200317'} 
        %datasetList = {'E20200317','N20180221'} 
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            case 'E20200317'
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',-50},{'state','Step 2','first',0}};
            case 'N20180221'
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',-50},{'state','Target Hold','first',0}};
        end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

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
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);
    
     
        
        minNumTimestamps = 10;
        
        %All Neuron PSTHs
        numGroups = ceil(numChannels/48);
        % Plot PSTH's
        for group = 1:numGroups
            neuronList = (group-1)*48+1:(group)*48;
            neuronList = neuronList(neuronList<=numChannels);
            f = figure; f.Position = [20 20 800 500];
            for posture = postureList
               for target = [1]
                  traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
                  time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
                  for neuron = neuronList
                      subplot(6,8,neuron-48*(group-1))
                      plot(time,traj(:,neuron),'-','Color',pcmap(posture,:));
                      hold on
                      if posture == postureList(end)
                        title(num2str(neuron))
                      end
                  end
               end
            end
            if saveFig
    %             saveas(gcf,fullfile(saveDir,[dataset,'_BCIPopTraj.svg']));
    %             saveas(gcf,fullfile(saveDir,[dataset,'_BCIPopTraj.fig']));
            end
        end
        
        %Select Neurons
        neuronList = [7,41,8,62,67,76];
        neuronList = [7,41,8];
        
        %Get min Time
        minTime = 10000;
        minNumTimestamps = 10000;
        for target = 1
            for posture = postureList
                time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
                if time(end) < minTime
                    minTime = time(end);
                    minNumTimestamps = length(time);
                end
            end
        end
        
        
        %Plot select PSTH's
        fs = 10;
        f = figure; f.Position = [10 10 350 500];
        [ha, pos] = tight_subplot(3,1,0.05,0.1,.15);
        maxTime = 0;
        numNeurons = length(neuronList);
        neuronInd = 1;
        for neuron = neuronList
           axes(ha(neuronInd));
           %Populate Figure
           for target = 1
               for posture = postureList
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
                    time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
                    if time(end) > maxTime
                        maxTime = time(end)
                    end
                    hold on
                    plot(time(1:minNumTimestamps),traj(1:minNumTimestamps,neuron),'Color',pcmap(posture,:),'LineWidth',2);
                    if posture == 5
                        title(['Monkey E Ch ',num2str(neuron)])
                    end
               end
           end
           if neuronInd == 2
                ylabel('FR (std)')
           end
           neuronInd = neuronInd + 1;
           ax = gca;
           ax.TickDir = 'out';
           xticks([0 250])
           ylimits = ax.YLim;
           yticks([ylimits(1) ylimits(2)])
           yticklabels({num2str(round(ylimits(1),1)), num2str(round(ylimits(2),1))})
           set(gca,'fontname','arial')
           set(gca,'fontsize',fs)
        end
        xlabel('time since go cue (ms)')
        xticklabels({'0','250'})
        for neuronInd = 1:3
            axes(ha(neuronInd));
            ax = gca;
            ax.XLim = [-50 minTime];
        end
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_BCI_PSTH.svg']));
            saveas(gcf,fullfile(saveDir,[dataset,'_BCI_PSTH.fig']));
        end
        
        
        %3D state-space view 
        posture = 5;
        target = 1;
        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj; 
        X = traj(:,neuronList(1:3));
 
        f = figure; %f.Position = [20 20 500 500];
        plot3(X(:,1),X(:,2),X(:,3),'-','Color',pcmap(posture,:),'LineWidth',2);
        hold on;
        plot3(X(1,1),X(1,2),X(1,3),'.','Color',pcmap(posture,:),'MarkerSize',20);
        plot3(X(end,1),X(end,2),X(end,3),'<','MarkerFaceColor',pcmap(posture,:),'MarkerEdgeColor',pcmap(posture,:),'MarkerSize',7);
        xlabel('Neuron 1 FR (std)')
        ylabel('Neuron 2 FR (std)')
        zlabel('Neuron 3 FR (std)')
        grid on
        axis equal
        ax = gca;
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_BCI_StateSpace.svg']));
            saveas(gcf,fullfile(saveDir,[dataset,'_BCI_StateSpace.fig']));
        end
        
        
                
        %Fit PC Plane 
        [coeff,score,latent,tsquared,explained,mu] = pca(X);
        planeColor = [.8,.8,.8];
        d = 0.25;
        

        
        %3D state-space view with Plane
        xlim = ax.XLim; ylim = ax.YLim;
        [surfX,surfY] = meshgrid(xlim(1):d:xlim(2),ylim(1):d:ylim(2));
        surfZ = nan(size(surfX));
        n = cross(coeff(:,1),coeff(:,2));
        r0 = mu;
        for row = 1:size(surfX,1)
            for col = 1:size(surfX,2)
                zX = surfX(row,col);
                zY = surfY(row,col);
                surfZ(row,col) = (dot(n,r0)-n(1)*zX-n(2)*zY)/n(3);
            end
        end
        %Z = x-mu*coeff;
        f = figure; %f.Position = [20 20 500 500];
        plot3(X(:,1),X(:,2),X(:,3),'-','Color',pcmap(posture,:),'LineWidth',2);
        hold on;
        plot3(X(1,1),X(1,2),X(1,3),'.','Color',pcmap(posture,:),'MarkerSize',20);
        plot3(X(end,1),X(end,2),X(end,3),'<','MarkerFaceColor',pcmap(posture,:),'MarkerEdgeColor',pcmap(posture,:),'MarkerSize',7);
        surf(surfX,surfY,surfZ,'EdgeColor','none','FaceColor',planeColor,'FaceAlpha',0.5);
        xlabel('Neuron 1 FR (std)')
        ylabel('Neuron 2 FR (std)')
        zlabel('Neuron 3 FR (std)')
        grid on
        axis equal
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_BCI_StateSpacePlane.svg']));
            saveas(gcf,fullfile(saveDir,[dataset,'_BCI_StateSpacePlane.fig']));
        end
        
        %Projection into plane 
        figure
        pcProj = (X-mu)*coeff;
        plot3(pcProj(:,1),pcProj(:,2),zeros(size(pcProj(:,1))),'-','Color',pcmap(posture,:),'LineWidth',2)
        hold on;
        plot3(pcProj(1,1),pcProj(1,2),0,'.','Color',pcmap(posture,:),'MarkerSize',20);
        plot3(pcProj(end,1),pcProj(end,2),0,'<','MarkerFaceColor',pcmap(posture,:),'MarkerEdgeColor',pcmap(posture,:),'MarkerSize',7);
        axis equal
        ax = gca;
        xlim = ax.XLim; ylim = ax.YLim;
        [surfX,surfY] = meshgrid([xlim(1),xlim(2)],[ylim(1),ylim(2)]);
        surfZ = zeros(size(surfX));
        surf(surfX,surfY,surfZ,'EdgeColor','none','FaceColor',planeColor,'FaceAlpha',0.5);
                plot3(pcProj(:,1),pcProj(:,2),zeros(size(pcProj(:,1))),'-','Color',pcmap(posture,:),'LineWidth',2)
        hold on;
        plot3(pcProj(1,1),pcProj(1,2),0,'.','Color',pcmap(posture,:),'MarkerSize',20);
        plot3(pcProj(end,1),pcProj(end,2),0,'<','MarkerFaceColor',pcmap(posture,:),'MarkerEdgeColor',pcmap(posture,:),'MarkerSize',7);
        
        xlabel(['PC 1 (',num2str(round(explained(1))),'%)'])
        ylabel(['PC 2 (',num2str(round(explained(2))),'%)'])
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        view([0 90])
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_BCI_StateSpaceProj.svg']));
            saveas(gcf,fullfile(saveDir,[dataset,'_BCI_StateSpaceProj.fig']));
        end
    end
    
    
    figure
    
