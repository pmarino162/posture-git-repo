clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\Interactions Reflect Task Demands';
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
                task = 'BCI';
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            case 'E20210706'
                task = 'Reach';
                %All
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',500}};
                %Reach
                reachTrialInclStates(1).trialName = {'GridReaching'};
                reachCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                reachTrialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',300}};
                %Center Hold
                centerHoldTrialInclStates(1).trialName = {'GridReaching'};
                centerHoldCondFields = {{'posture','conditionData','postureID'}};
                centerHoldTrialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Delay','first',50}};
                %Use Reach
                trialInclStates = reachTrialInclStates;
                condFields = reachCondFields;
                
                
        end
        
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        if strcmp(task,'Reach')
           trajStruct = trajStruct(ismember([trajStruct.posture],[1:4]));
           centerHoldTrajStruct = getTrajStruct20220419(Data,centerHoldCondFields,trajFields,centerHoldTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
           centerHoldTrajStruct = centerHoldTrajStruct(ismember([centerHoldTrajStruct.posture],[1:4]));
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
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);
    

        %Fit Posture LDA
        obsStruct = struct('label',[],'numObs',[],'allObs',[]);
        structInd = 1;   
        for posture = postureList
           obsStruct(structInd).label = posture;
           allObs = [];
           if strcmpi(task,'BCI')
               tempTrajStruct = trajStruct([trajStruct.posture]==posture);
           elseif strcmpi(task,'Reach')
               tempTrajStruct = centerHoldTrajStruct([centerHoldTrajStruct.posture]==posture);
           end
           for i = 1:size(tempTrajStruct,2)
               for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
                   timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
                   traj = tempTrajStruct(i).allSmoothFR(j).traj;
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
        [postureLDA] = doLDA(obsStruct);
        
        %Fit BCI PC
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
        allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
        [allPCA,score,latent,tsquared,explained,allMu] = pca(allTraj);
        
        %Orthonormalize BCI Space 
        if strcmpi(task,'Reach')
            [orthProj,~] = qr([postureLDA(:,1),allPCA(:,1),allPCA(:,3)]); orthProj = orthProj(:,1:3);
        elseif strcmpi(task,'BCI')
            [orthProj,~] = qr([postureLDA(:,1),allPCA(:,1),allPCA(:,3)]); orthProj = orthProj(:,1:3);
        end
        
        %Make sure that posture 1 is bottom and posture 4 is top
        P1Proj = mean(obsStruct([obsStruct.label]==1).allObs)*orthProj;
        P4Proj = mean(obsStruct([obsStruct.label]==3).allObs)*orthProj;
        if P1Proj(:,1) > P4Proj(:,1)
            orthProj(:,1) = -1.*orthProj(:,1);
        end
        
        
        %Add Projections to trajStruct
        totalVar = trace(cov(allTraj));
        for i = 1:size(trajStruct,2)
            %Add average traces to trajStruct
            trajStruct(i).orthProj.traj = trajStruct(i).avgSmoothFR.traj*orthProj;
            %Get VAF
            trajStruct(i).orthProj.VAF =  100.*(diag(cov(allTraj*orthProj))')./totalVar;
        end
        VAF = round(trajStruct(1).orthProj.VAF); 
        
        %Plot groups of trajectories for each task
        %Plot full data and get axis limits/viewing angles
        fs = 14;
        condGroup = struct('postureList',[],'targetList',[]);
        condGroup(1).postureList = [1:4];
        if strcmp(task,'Reach')
            condGroup(1).targetList = 2;
        elseif strcmp(task,'BCI')
            condGroup(1).targetList = 1;
        end
        for condGroupInd = 1
            postureList = condGroup(condGroupInd).postureList;
            targetList = condGroup(condGroupInd).targetList;
            f = figure; set(gcf,'color','white'); hold on
            hold on
            timePts = 1:minNumTimestamps;
            xDim = 2; yDim = 3; zDim = 1;
            for posture = postureList
                for target = targetList
                    if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).orthProj.traj(timePts,:); 
                        plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    end
                end
            end
            
            
            grid on
            xticklabels({}); yticklabels({}); zticklabels({}); 
            %xlabel(['Neural Dim 1 (',num2str(VAF(xDim)),'%)']); ylabel(['Neural Dim 2 (',num2str(VAF(yDim)),'%)']); zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
            xlabel(['Neural Dim 1']); ylabel(['Neural Dim 2']); zlabel(['Posture Dim 1'])
            view([20 10])
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
            %axis equal
            axis vis3d
            saveas(gcf,fullfile(saveDir,[task,'.svg']));
        end
        
        ax = gca; xlimits = ax.XLim; ylimits = ax.YLim; zlimits = ax.ZLim;
        
        %Save top-down view highlighting trajectory differences
        view([90 90])
        saveas(gcf,fullfile(saveDir,[task,'TrajDifferences.svg']));
        view([20 10])
        
        %Animation 1 = Rotate full plot 
        fileName = fullfile(saveDir,[dataset,'rotating.mp4']);
        xlabel([]); ylabel([]); zlabel([])
        
        OptionZ.FrameRate=12;OptionZ.Duration=5.5;OptionZ.Periodic=true;
        %viewVec = [90,10; 0,10; -90,10; -180,10; -270,10];
        %viewVec = [20,10; -70,10; -160,10; -250,10; -340,10];
        viewVec = [20,10; 0,90];
        CaptureFigVid(viewVec,fileName,OptionZ)
        

        
        
        %Animation 2 = Keep previous trajectory; animate new trajectory
        fps = 14;
        f = figure; set(gcf,'color','white'); hold on
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1;
        postureList = condGroup(1).postureList; targetList = condGroup(1).targetList;
        for posture = postureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    fileName = fullfile(saveDir,[dataset,'trajSequence_','P',num2str(posture),'.mp4']);
                    v = VideoWriter(fileName); v.FrameRate = fps; open(v);
                    frame = 1;
                    if frame == 1;
                        plot3(0,0,0);
                        drawnow
                        xlim(xlimits); ylim(ylimits); zlim(zlimits); grid on; xticklabels({}); yticklabels({}); zticklabels({}); 
                        %xlabel(['Neural Dim 1 (',num2str(VAF(xDim)),'%)']); ylabel(['Neural Dim 2 (',num2str(VAF(yDim)),'%)']); zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
                        xlabel(['Neural Dim 1']); ylabel(['Neural Dim 2']); zlabel(['Posture Dim 1'])
                        view([20 10])
                        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
                        M(frame) = getframe(gcf);writeVideo(v,M(frame));
                        frame = frame + 1;
                    end
                    for t = timePts
                        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).orthProj.traj(timePts,:); 
                        plot3(traj(1:t,xDim),traj(1:t,yDim),traj(1:t,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                        if t == timePts(end)
                            plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                        end
                        drawnow
                        xlim(xlimits); ylim(ylimits); zlim(zlimits); grid on; xticklabels({}); yticklabels({}); zticklabels({}); 
                        %xlabel(['Neural Dim 1 (',num2str(VAF(xDim)),'%)']); ylabel(['Neural Dim 2 (',num2str(VAF(yDim)),'%)']); zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
                        xlabel(['Neural Dim 1']); ylabel(['Neural Dim 2']); zlabel(['Posture Dim 1'])
                        view([20 10])
                        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
                        M(frame) = getframe(gcf); writeVideo(v,M(frame));
                        frame = frame + 1;
                    end
                    close(v);
                end
            end
        end
        
        
        
        
        %Animation 3 = Animate trajectories in sequence
        fps = 14;
        fileName = fullfile(saveDir,[dataset,'trajSequence.mp4']);
        v = VideoWriter(fileName); v.FrameRate = fps; open(v);
        postureList = condGroup(1).postureList; targetList = condGroup(1).targetList;
        f = figure; set(gcf,'color','white'); hold on
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1;
        frame = 1;
        for posture = postureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    for t = timePts
                        if frame ~= 1
                            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).orthProj.traj(timePts,:); 
                            plot3(traj(1:t,xDim),traj(1:t,yDim),traj(1:t,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                            if t == timePts(end)
                                plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                            end
                        end
                        drawnow
                        xlim(xlimits); ylim(ylimits); zlim(zlimits); grid on; xticklabels({}); yticklabels({}); zticklabels({}); 
                        %xlabel(['Neural Dim 1 (',num2str(VAF(xDim)),'%)']); ylabel(['Neural Dim 2 (',num2str(VAF(yDim)),'%)']); zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
                        xlabel(['Neural Dim 1']); ylabel(['Neural Dim 2']); zlabel(['Posture Dim 1'])
                        view([20 10])
                        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
                        M(frame) = getframe(gcf); writeVideo(v,M(frame));
                        frame = frame + 1;
                    end
                end
            end
        end
        close(v);
        movie(M,1,fps)
        
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