clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
     saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\BCIPopTraj\Animations';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
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
        binWidth = 50;
        kernelStdDev = 50;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            case 'E20200317'
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            case 'N20180221'
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
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
    
        %Get pcmap for number of postures
%         pcmap = orli(numPostures);
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
    
        %Do PCA on all condition averages
        allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
        [allPCA,score,latent,tsquared,explained,allMu] = pca(allTraj);

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
           tempTrajStruct = trajStruct([trajStruct.posture]==posture);
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
        
        %Peform LDA by target on all data 
        obsStruct = struct('label',[],'numObs',[],'allObs',[]);
        structInd = 1;   
        for target = targetList
           obsStruct(structInd).label = target;
           allObs = [];
           tempTrajStruct = trajStruct([trajStruct.target]==target);
           timeToUse = 125;
           numOnEitherSide = 1;
           for i = 1:size(tempTrajStruct,2)
               for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
                   timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
                   traj = tempTrajStruct(i).allSmoothFR(j).traj;
                   %[~,timeToUseInd] = min(abs(timestamps-timeToUse));
                   timeToUseInd = size(traj,1)-1;
                   traj = traj(timeToUseInd-numOnEitherSide:timeToUseInd+numOnEitherSide,:);
                   allObs = vertcat(allObs,traj);
%                    if size(traj,1) > minNumTimestamps
%                         traj = traj(1:minNumTimestamps,:);
%                    end
                   allObs = vertcat(allObs,traj);
               end
           end
           obsStruct(structInd).allObs = allObs;
           obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
           structInd = structInd + 1;
        end
        [targetLDA] = doLDA(obsStruct);
        

        
        %Orthonormalize combinations of axes
        [CPTOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1),CIPCA(:,1)]); CPTOrth = CPTOrth(:,1:3);
        [CPOrth,~] = qr([postureLDA(:,1:2),CIPCA(:,1)]); CPOrth = CPOrth(:,1:3);
        [PTOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
        [PTOrthLDA,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]); PTOrthLDA = PTOrthLDA(:,1:3);
        [PTOrthPCA,~] = qr([posturePCA(:,1),targetNoCPPCA(:,1:2)]); PTOrthPCA = PTOrthPCA(:,1:3);
        [CPTAllOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1:2),CIPCA(:,1)]); CPTAllOrth = CPTAllOrth(:,1:4);

        %Flip Posture Dim so that posture 1 is on the left
        p1TrajMean = mean(trajStruct([trajStruct.posture]==1 & [trajStruct.target]==1).avgSmoothFR.traj,1); 
        p2TrajMean = mean(trajStruct([trajStruct.posture]==max(postureList) & [trajStruct.target]==1).avgSmoothFR.traj,1); 
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
        
        
    
        %Add Projections to trajStruct
        totalVar = trace(cov(allTraj));
        for i = 1:size(trajStruct,2)
            %Add average traces to trajStruct
            trajStruct(i).PCA.traj = (trajStruct(i).avgSmoothFR.traj-allMu)*allPCA;
            trajStruct(i).CIPCA.traj = (trajStruct(i).avgSmoothFR.traj-CIMargOffset')*CIPCA;
            trajStruct(i).targetPCA.traj = trajStruct(i).avgSmoothFR.traj*targetPCA;
            trajStruct(i).posturePCA.traj = trajStruct(i).avgSmoothFR.traj*posturePCA;
            trajStruct(i).postureLDA.traj = trajStruct(i).avgSmoothFR.traj*postureLDA;
            trajStruct(i).CPTOrth.traj = trajStruct(i).avgSmoothFR.traj*CPTOrth;
            trajStruct(i).CPOrth.traj = trajStruct(i).avgSmoothFR.traj*CPOrth;
            trajStruct(i).PTOrth.traj = trajStruct(i).avgSmoothFR.traj*PTOrth;
            trajStruct(i).PTOrthLDA.traj = trajStruct(i).avgSmoothFR.traj*PTOrthLDA;
            trajStruct(i).PTOrthPCA.traj = trajStruct(i).avgSmoothFR.traj*PTOrthPCA;
            trajStruct(i).CPTAllOrth.traj = trajStruct(i).avgSmoothFR.traj*CPTAllOrth;
            %Get VAF
            trajStruct(i).PCA.VAF =  100.*(diag(cov(allTraj*allPCA))')./totalVar;
            trajStruct(i).PTOrth.VAF =  100.*(diag(cov(allTraj*PTOrth))')./totalVar;
            trajStruct(i).PTOrthLDA.VAF =  100.*(diag(cov(allTraj*PTOrthLDA))')./totalVar;
            trajStruct(i).PTOrthPCA.VAF =  100.*(diag(cov(allTraj*PTOrthPCA))')./totalVar;
            trajStruct(i).postureLDA.VAF = 100.*(diag(cov(allTraj*postureLDA))')./totalVar;
        end
        
        minNumTimestamps = 6;
        fs = 14;
        
        %Plot full data and get axis limits/viewing angles
        condGroup = struct('postureList',[],'targetList',[]);
        condGroup(1).postureList = [1:5];
            condGroup(1).targetList = 1:8;
        for condGroupInd = 1
            postureList = condGroup(condGroupInd).postureList;
            targetList = condGroup(condGroupInd).targetList;
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
            
            if condGroupInd == 1
                ax = gca; 
                xlimits = ax.XLim;
                ylimits = ax.YLim;
                zlimits = ax.ZLim;
            end
            xlim(xlimits)
            ylim(ylimits)
            zlim(zlimits)
            
            grid on
            xticklabels({}); yticklabels({}); zticklabels({}); 
            VAF = round(trajStruct(1).PTOrth.VAF);

            xlabel(['Target Dim 1 (',num2str(VAF(xDim)),'%)'])
            ylabel(['Target Dim 2 (',num2str(VAF(yDim)),'%)']) 
            zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
            view([20 10])
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
        end
        
        
        %Set animation params 
        fps = 10;
        
        %Plot animation of P1T1 proceeding over time
        fileName = fullfile(saveDir,[dataset,'_P1T1.mp4']);
        v = VideoWriter(fileName); v.FrameRate = fps; open(v);
        postureList = 1; targetList = 1;
        f = figure; set(gcf,'color','white'); hold on
        timePts = 1:minNumTimestamps;
        frameList = [1,timePts+1];
        xDim = 2; yDim = 3; zDim = 1; VAF = round(trajStruct(1).PTOrth.VAF);
        for frame = frameList
            for posture = postureList
                for target = targetList
                    if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                        if frame ~= 1
                            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                            plot3(traj(1:frame-1,xDim),traj(1:frame-1,yDim),traj(1:frame-1,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                            if frame == frameList(end)
                                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                            end
                        end
                        drawnow
                    end
                end
            end
            xlim(xlimits); ylim(ylimits); zlim(zlimits); grid on; xticklabels({}); yticklabels({}); zticklabels({}); 
            xlabel(['Target Dim 1 (',num2str(VAF(xDim)),'%)']); ylabel(['Target Dim 2 (',num2str(VAF(yDim)),'%)']); zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
            view([20 10])
            set(gca,'fontname','arial'); set(gca,'fontsize',fs)
            M(frame) = getframe(gcf); writeVideo(v,M(frame));
        end
        close(v);
        movie(M,1,fps)
        
        %Keep P1T1; Add rest of targets
        fileName = fullfile(saveDir,[dataset,'_P1allT.mp4']);
        v = VideoWriter(fileName); v.FrameRate = fps; open(v);
        condGroup = struct('postureList',[],'targetList',[]);
        f = figure; set(gcf,'color','white'); hold on
        postureList = 1; targetList = 1;
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1; VAF = round(trajStruct(1).PTOrth.VAF);
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
        frameList = [1,timePts+1];
        postureList = 1; targetList = 2:8;     
        for frame = frameList
            for posture = postureList
                for target = targetList
                    if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                        if frame ~= 1
                            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                            plot3(traj(1:frame-1,xDim),traj(1:frame-1,yDim),traj(1:frame-1,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                            if frame == frameList(end)
                                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                            end
                        end
                        drawnow 
                    end
                end
            end
            xlim(xlimits); ylim(ylimits); zlim(zlimits); grid on; xticklabels({}); yticklabels({}); zticklabels({}); 
            xlabel(['Target Dim 1 (',num2str(VAF(xDim)),'%)']); ylabel(['Target Dim 2 (',num2str(VAF(yDim)),'%)']); zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
            view([20 10])
            set(gca,'fontname','arial'); set(gca,'fontsize',fs)
             M(frame) = getframe(gcf); writeVideo(v,M(frame));
        end
        close(v);
        movie(M,1,fps)
        
        
        %Keep P1; Add P2
        fileName = fullfile(saveDir,[dataset,'_P1andP2.mp4']);
        v = VideoWriter(fileName); v.FrameRate = fps; open(v);
        condGroup = struct('postureList',[],'targetList',[]);
        f = figure; set(gcf,'color','white'); hold on
        postureList = 1; targetList = 1:8;
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1; VAF = round(trajStruct(1).PTOrth.VAF);
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
        frameList = [1,timePts+1];
        postureList = 2; targetList = 1:8;     
        for frame = frameList
            for posture = postureList
                for target = targetList
                    if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                        if frame ~= 1
                            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                            plot3(traj(1:frame-1,xDim),traj(1:frame-1,yDim),traj(1:frame-1,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                            if frame == frameList(end)
                                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                            end
                        end
                        drawnow 
                    end
                end
            end
            xlim(xlimits); ylim(ylimits); zlim(zlimits); grid on; xticklabels({}); yticklabels({}); zticklabels({}); 
            xlabel(['Target Dim 1 (',num2str(VAF(xDim)),'%)']); ylabel(['Target Dim 2 (',num2str(VAF(yDim)),'%)']); zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
            view([20 10])
            set(gca,'fontname','arial'); set(gca,'fontsize',fs)
             M(frame) = getframe(gcf); writeVideo(v,M(frame));
        end
        close(v);
        movie(M,1,fps)
        
        
        
        %Keep P1 and P2; add other postures
        fileName = fullfile(saveDir,[dataset,'_restOfPostures.mp4']);
        v = VideoWriter(fileName); v.FrameRate = fps; open(v);
        condGroup = struct('postureList',[],'targetList',[]);
        f = figure; set(gcf,'color','white'); hold on
        postureList = 1:2; targetList = 1:8;
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1; VAF = round(trajStruct(1).PTOrth.VAF);
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
        frameList = [1,timePts+1];
        postureList = 3:5; targetList = 1:8;     
        for frame = frameList
            for posture = postureList
                for target = targetList
                    if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                        if frame ~= 1
                            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                            plot3(traj(1:frame-1,xDim),traj(1:frame-1,yDim),traj(1:frame-1,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                            if frame == frameList(end)
                                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                            end
                        end
                        drawnow 
                    end
                end
            end
            xlim(xlimits); ylim(ylimits); zlim(zlimits); grid on; xticklabels({}); yticklabels({}); zticklabels({}); 
            xlabel(['Target Dim 1 (',num2str(VAF(xDim)),'%)']); ylabel(['Target Dim 2 (',num2str(VAF(yDim)),'%)']); zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
            view([20 10])
            set(gca,'fontname','arial'); set(gca,'fontsize',fs)
             M(frame) = getframe(gcf); writeVideo(v,M(frame));
        end
        close(v);
        movie(M,1,fps)
        
        
        %Rotate full plot 
        fileName = fullfile(saveDir,[dataset,'_RotatingSpiders.mp4']);
        condGroup = struct('postureList',[],'targetList',[]);
        f = figure; set(gcf,'color','white'); hold on
        postureList = 1:5; targetList = 1:8;
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1; VAF = round(trajStruct(1).PTOrth.VAF);
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
        xlim(xlimits); ylim(ylimits); zlim(zlimits); grid on; xticklabels({}); yticklabels({}); zticklabels({}); 
        %xlabel(['Target Dim 1']); ylabel(['Target Dim 2']); zlabel(['Posture Dim 1']);
        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        p = get(gca,'Position');
        set(gca,'Position',[p(1) p(2)+0.1 p(3)-0.1 p(4)-0.2]);
        axis vis3d
        OptionZ.FrameRate=12;OptionZ.Duration=5.5;OptionZ.Periodic=true;
        %viewVec = [90,10; 0,10; -90,10; -180,10; -270,10];
        viewVec = [20,10; -70,10; -160,10; -250,10; -340,10];
        CaptureFigVid(viewVec,fileName,OptionZ)
    
        
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
         