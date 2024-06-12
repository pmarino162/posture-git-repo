clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 2';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set parameters
    regularize = true;
    
%% Main Loop  
    for datasetList = {'E20210708'}%{'E20200318','E20210901','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}
        %% Load Data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);

        %% Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);
    
        %% For earl reaching, keep postures 2 and 7, but relabel as 1 and 2 for color purposes
        if strcmpi(dataset,'E20210708')
            trajStruct = trajStruct(ismember([trajStruct.posture],[2,7]));
            for i = 1:size(trajStruct,2)
               if trajStruct(i).posture == 2
                   trajStruct(i).posture = 1;
               elseif trajStruct(i).posture == 7
                   trajStruct(i).posture = 2;
               end
            end
        end
              
        %% Get traj struct dimensions
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 
        
        %% Setup colormap (based on number of postures)
        [pcmap,tcmap,rainbow] = getColorMaps(numPostures);  

        %% Run dPCA
        [W,V,postureDims,targetDims,explVar,additionalVarExpl] = getPostureDPCADims(trajStruct,regularize,postureList,numPostures,targetList,numTargets,numChannels,minNumTimestamps);
        pDPCA = W(:,postureDims);
        tDPCA = W(:,targetDims);     
        [PTOrth,~] = qr([pDPCA(:,1),tDPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
                     
        %Add Projections to trajStruct
        for i = 1:size(trajStruct,2)
            %Add average traces to trajStruct
            trajStruct(i).PTOrth.traj = trajStruct(i).avgZSmoothFR.traj*PTOrth;
            trajStruct(i).dPCA.traj = trajStruct(i).avgZSmoothFR.traj*W;
            %Get VAF
            trajStruct(i).PTOrth.VAF =  [additionalVarExpl(postureDims(1)),additionalVarExpl(targetDims(1)),additionalVarExpl(targetDims(2))];
            trajStruct(i).dPCA.VAF = explVar.componentVar;
        end

        %Plot - Orthographic
        fs = 14;
        figure
        hold on
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1;
        for posture = postureList
            for target = targetList
                 if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                    plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
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
        end
        
        %Plot vs time
        figure 
        timePts = 1:minNumTimestamps;
        for posture = postureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).dPCA.traj(timePts,:); 
                    time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgZSmoothFR.timestamps(timePts); 
                    for dim = 1:10
                       subplot(2,5,dim)
                       plot(time,traj(:,dim),'Color',pcmap(posture,:),'LineWidth',2);
                       hold on;                        
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
            VAF = round(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).dPCA.VAF);
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

  %      clf; close all
    end