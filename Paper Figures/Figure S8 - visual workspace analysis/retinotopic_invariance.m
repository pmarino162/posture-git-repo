clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = "C:\Users\pmari\OneDrive\Documents\Posture\Paper\Reviewer responses\Analyses\Retinotopic invariance";
    set(0, 'DefaultFigureRenderer', 'painters');

%% Set parameters
    regularize = true;

%% Load N and R data (all visual data)
    %reachDatasetList = {'N20190222_all_visual_data','N20190226_all_visual_data',...
    %    'N20190227_all_visual_data','N20190228_all_visual_data','N20190307_all_visual_data'...
    %    'R20200221_all_visual_data','R20200222_all_visual_data'}; 

    reachDatasetList = {'R20200221_all_visual_data'};
    
%% Main loop   
    for datasetList = reachDatasetList
        % Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);                 
        % Split trajStruct based on visualID
        visualID = [trajStruct.visual];
        trajStructV1 = trajStruct(visualID == 1);
        trajStructV2 = trajStruct(visualID == 2);
        % Get trajStruct dimensions for visualID = 2
        [minNumTimestamps] = getMinNumTimestamps(trajStructV2); 
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStructV2); 
        % Run dPCA on visualID = 2
        [W,V,postureDims,targetDims,explVar,additionalVarExpl] = getPostureDPCADims(trajStructV2,regularize,postureList,numPostures,targetList,numTargets,numChannels,minNumTimestamps);
        pDPCA = W(:,postureDims);
        tDPCA = W(:,targetDims);     
        [PTOrth,~] = qr([pDPCA(:,1),tDPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
        
        % Add projections for both visualIDs to trajStruct
        for i = 1:size(trajStruct,2)
            %Add average traces to trajStruct
            trajStruct(i).PTOrth.traj = trajStruct(i).avgZSmoothFR.traj*PTOrth;
            trajStruct(i).dPCA.traj = trajStruct(i).avgZSmoothFR.traj*W;
            %Get VAF
            trajStruct(i).PTOrth.VAF =  [additionalVarExpl(postureDims(1)),additionalVarExpl(targetDims(1)),additionalVarExpl(targetDims(2))];
            trajStruct(i).dPCA.VAF = explVar.componentVar;
        end
        
        
        
        
    end
    
    
%% Plot
        %% Setup colormap (based on number of postures)
        [pcmap,tcmap,rainbow] = getColorMaps(numPostures);  

        %Plot - Orthographic
        fs = 14;
        figure
        hold on
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1;
        for visualID = [1,2]
            for posture = postureList
                for target = targetList
                     if any([trajStruct.visual]==visualID & [trajStruct.posture]==posture & [trajStruct.target]==target)
                        traj = trajStruct([trajStruct.visual]==visualID & [trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                        if visualID == 1
                            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                        elseif visualID == 2
                            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'--','Color',pcmap(posture,:),'LineWidth',2);
                        end
                        
                        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    end
                end
            end
        end
             
        grid on
        xticklabels({}); yticklabels({}); zticklabels({}); 
        VAF = round(trajStruct(1).PTOrth.VAF);

        xlabel(['Goal Dim 1 (',num2str(VAF(xDim)),'%)'])
        ylabel(['Goal Dim 2 (',num2str(VAF(yDim)),'%)']) 
        zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
        view([20 10])
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PTOrth.fig']));
            saveas(gcf,fullfile(saveDir,[dataset,'_PTOrth.svg']));
        end
        
        
        % Plot - 2d goal space
        fs = 14;

        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1;
        

        for posture = postureList
            figure
            hold on
            for visualID = [1,2]
                for target = targetList
                     if any([trajStruct.visual]==visualID & [trajStruct.posture]==posture & [trajStruct.target]==target)
                        traj = trajStruct([trajStruct.visual]==visualID & [trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                        if visualID == 1
                            plot(traj(:,xDim),traj(:,yDim),'Color',pcmap(posture,:),'LineWidth',2);
                        elseif visualID == 2
                            plot(traj(:,xDim),traj(:,yDim),'--','Color',pcmap(posture,:),'LineWidth',2);
                        end
                        
                        plot(traj(1,xDim),traj(1,yDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                        plot(traj(end,xDim),traj(end,yDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    end
                end
            end
            
            grid on
            xticklabels({}); yticklabels({}); zticklabels({}); 
            VAF = round(trajStruct(1).PTOrth.VAF);

            xlabel(['Goal Dim 1 (',num2str(VAF(xDim)),'%)'])
            ylabel(['Goal Dim 2 (',num2str(VAF(yDim)),'%)']) 
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_posture_',num2str(posture),'_GoalSubspace.jpg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_posture_',num2str(posture),'_GoalSubspace.svg']));
            end
            
        end
             

        
        
        %Plot vs time
        figure 
        timePts = 1:minNumTimestamps;
        for posture = postureList
            for target = [1,5]
                for visualID = [1,2]
                    if any([trajStruct.visual]==visualID & [trajStruct.posture]==posture & [trajStruct.target]==target)
                        traj = trajStruct([trajStruct.visual]==visualID & [trajStruct.posture]==posture & [trajStruct.target]==target).dPCA.traj(timePts,:); 
                        time = trajStruct([trajStruct.visual]==visualID & [trajStruct.posture]==posture & [trajStruct.target]==target).avgZSmoothFR.timestamps(timePts); 
                        for dim = 1:10
                           subplot(2,5,dim)
                           if visualID == 1
                               plot(time,traj(:,dim),'Color',pcmap(posture,:),'LineWidth',2);
                           elseif visualID == 2
                               plot(time,traj(:,dim),'--','Color',pcmap(posture,:),'LineWidth',2);
                           end
                           hold on;                        
                        end
                    end
                end
            end
        end
          
        %Format
        maxYRange = 0;
        for dim = 1:10
            subplot(2,5,dim)
            ax = gca;
            ylimits = ax.YLim;
            yRange = ylimits(2)-ylimits(1);
            if yRange > maxYRange
                maxYRange = yRange;
            end
        end
        for dim = 1:10
            VAF = round(trajStruct([trajStruct.visual]==visualID & [trajStruct.posture]==posture & [trajStruct.target]==target).dPCA.VAF);
            subplot(2,5,dim);
            ax = gca;
            ylimits = ax.YLim;
            yMid = ylimits(1) + 0.5*(ylimits(2)-ylimits(1));
            ylim([yMid-maxYRange/2 yMid+maxYRange/2])
            ylabel(['dPC ',num2str(dim),' (',num2str(VAF(dim)),'%)'])
            yticks([min(ax.YTick) max(ax.YTick)]);
            ax.TickDir = 'out';
            if dim < 4
                xticks([])
            else
                xticks([-200 0 300]);
            end
            if dim == 5
                %xlabel('time (ms)')
            end
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
        end    