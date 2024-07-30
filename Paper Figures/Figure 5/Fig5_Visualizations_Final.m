clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 5';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Parameters
    numPCsToKeep = 10; %Num PCs to project data into before analysis
    pctCI = 95;
    ellAlpha = 0.4;
    
%% Load and preprocess data
    %Load data
    dataset = 'E20200314';
    [Data,zScoreParams] = loadData(dataset);

    %Get trajStruct
    [Data] = removeShortBCIandIsoTrials(Data,dataset);
    [~,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);  
    [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 
    taskList = 1:3;
    % Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(numPostures);
    %For each timecourse, get one point
     for i = 1:size(trajStruct,2)
        for j = 1:size(trajStruct(i).allZSmoothFR,2)
           trajStruct(i).allZSmoothFR(j).obs = mean(trajStruct(i).allZSmoothFR(j).traj); 
        end
     end
    %Reduce dimensionality of observations
    allObs = [];
    for i = 1:size(trajStruct,2)
        for j = 1:size(trajStruct(i).allZSmoothFR,2)
            obs = trajStruct(i).allZSmoothFR(j).obs;
            allObs = vertcat(allObs,obs);
        end
    end
    [coeff,score,latent,tsquared,explained,mu] = pca(allObs); 
    for i = 1:size(trajStruct,2)
        for j = 1:size(trajStruct(i).allZSmoothFR,2)
            trajStruct(i).allZSmoothFR(j).obs = (trajStruct(i).allZSmoothFR(j).obs-mu)*coeff(:,1:numPCsToKeep);
        end
    end
    
        
 %% Create training set, doing appropriate balancing
    %Balancing scheme: for across-task posture decoding, balance by targets
    %w/in each posture x task, and get same number of points from each task. (this scheme doesn't balance perfectly, but it also keeps more data) 
    %For task-decoding, balance by target and posture
    %w/in each task
 
    %Get minimum number of trials to any target direction within each
    %task x posture. This can be used for both across-task posture decoding
    %and task decoding. 
    minNumObsStruct = struct('task',[],'posture',[],'minNumObs',[]);
    structInd = 1;
    for task = 1:3
        for posture = 1:3
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            tempMinNumCondObs = 999999999;
            for i = 1:size(tempTrajStruct,2)
                numCondObs = size(vertcat(tempTrajStruct(i).allZSmoothFR.obs),1);
                if numCondObs < tempMinNumCondObs
                    tempMinNumCondObs = numCondObs;
                end
            end
            minNumObsStruct(structInd).task = task;
            minNumObsStruct(structInd).posture = posture;
            minNumObsStruct(structInd).minNumObs = tempMinNumCondObs;
            structInd = structInd + 1;
        end
    end
    
    %Build taskObsStruct, which is balanced by target for each task
    taskObsStruct = struct('task',[],'obs',[],'labels',[]);
    structInd = 1;
    for task = 1:4
        obs = []; labels = [];
        if task < 4
            for posture = 1:3
                tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
                numObsToUsePerTarget = minNumObsStruct([minNumObsStruct.task]==task & [minNumObsStruct.posture]==posture).minNumObs;
                for i = 1:size(tempTrajStruct,2)
                    condObs = vertcat(tempTrajStruct(i).allZSmoothFR.obs);
                    obs = vertcat(obs,datasample(condObs,numObsToUsePerTarget,'Replace',false));
                    posture = tempTrajStruct(i).posture;
                    labels = vertcat(labels,ones(numObsToUsePerTarget,1)*posture);
                end  
            end
        elseif task == 4
           minNumTaskObs = 999999999999;
           for tempTask = 1:3
               numTaskObs = size(taskObsStruct(tempTask).obs,1);
               if numTaskObs < minNumTaskObs
                   minNumTaskObs = numTaskObs;
               end
           end
           for tempTask = 1:3
                taskObs = taskObsStruct(tempTask).obs;
                taskLabels = taskObsStruct(tempTask).labels;
                [tempObs,idx] = datasample(taskObs,minNumTaskObs,'Replace',false);
                tempLabels = taskLabels(idx);
                obs = vertcat(obs,tempObs);
                labels = vertcat(labels,tempLabels);
           end
        end
        taskObsStruct(structInd).task = task;
        taskObsStruct(structInd).obs = obs;
        taskObsStruct(structInd).labels = labels;
        structInd = structInd + 1; 
    end
    
    
    %Build training set 
    trainingSets = struct('decoder','','obs',[],'labels',[]);
    for decoder = 1:2
        obs = []; labels = [];
        if decoder == 1 %Across-task posture decoder
            trainingSets(decoder).decoder = 'posture';
            obs = taskObsStruct([taskObsStruct.task]==4).obs;
            labels = taskObsStruct([taskObsStruct.task]==4).labels;            
        elseif decoder == 2
            %For each task, find posture x target condition with lowest
            %trial count. Sample that number of trials for every condition
            %within that task.
            for task = 1:3
                minNumObs = min([minNumObsStruct([minNumObsStruct.task]==task).minNumObs]);
                for posture = 1:3
                    tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
                    for i = 1:size(tempTrajStruct,2) %Loop over targets
                        condObs = vertcat(tempTrajStruct(i).allZSmoothFR.obs);
                        tempObs = datasample(condObs,minNumObs,'Replace',false);
                        tempLabels = task*ones(minNumObs,1);
                        obs = vertcat(obs,tempObs);
                        labels = vertcat(labels,tempLabels);   
                    end
                end
            end
            trainingSets(decoder).decoder = 'task';
        end
        trainingSets(decoder).obs = obs;
        trainingSets(decoder).labels = labels;
    end


%% Do LDA by posture and task on all trial avgerages
    obs = trainingSets(1).obs; labels = trainingSets(1).labels;
    numClasses = 3;
    LDAProj = fisherLDA(obs, labels);
    [LDAProj,~] = qr(LDAProj);
    LDAProj = LDAProj(:,1:numClasses-1);
    postureLDA = LDAProj;
    
    %Flip axes so that P1 is on the right and P2 is on the top 
    P1Proj = mean(obs(labels==1,:)*postureLDA);
    P2Proj = mean(obs(labels==2,:)*postureLDA);
    if P1Proj(:,1) < P2Proj(:,1)
        postureLDA(:,1) = -1.*postureLDA(:,1);
    end
    if P1Proj(:,2) > P2Proj(:,2)
        postureLDA(:,2) = -1.*postureLDA(:,2);
    end
  
    
%% Do LDA by task
    obs = trainingSets(2).obs; labels = trainingSets(2).labels;
    numClasses = 3;
    LDAProj = fisherLDA(obs, labels);
    [LDAProj,~] = qr(LDAProj);
    LDAProj = LDAProj(:,1:numClasses-1);
    taskLDA = LDAProj;
    
    %Choose task dim within 2d task space
    taskDim = taskLDA(:,1)-0.5.*taskLDA(:,2);
    taskDim = taskDim./vecnorm(taskDim);
    
    [PTaskOrth,~] = qr([postureLDA,taskDim]); PTaskOrth = PTaskOrth(:,1:3);
    %Make sure posture dims match postrureLDA directions
    PTaskOrth = [postureLDA,PTaskOrth(:,3)];
%     %Flip axes so that P1 is on the right and P2 is on the top 
%     P1Proj = mean(obs(labels==1,:)*PTaskOrth);
%     P2Proj = mean(obs(labels==2,:)*PTaskOrth);
%     if P1Proj(:,1) < P2Proj(:,1)
%         PTaskOrth(:,1) = -1.*PTaskOrth(:,1);
%     end
%     if P1Proj(:,2) > P2Proj(:,2)
%         PTaskOrth(:,2) = -1.*PTaskOrth(:,2);
%     end
%     

%% Get posture error ellipses and means for all tasks
    errorEllipses = struct('task',[],'posture',[],'X',[],'Y',[],'u',[]);
    structInd = 1;
    for task = taskList
        for posture = postureList
            allTrialAvg = [];
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            for i = 1:size(tempTrajStruct,2)
                allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allZSmoothFR.obs);
            end
            ldaProj = allTrialAvg*postureLDA;
            [X,Y,Z] = error_ellipse_Patrick(ldaProj,pctCI);
            errorEllipses(structInd).task = task;
            errorEllipses(structInd).posture = posture;
            errorEllipses(structInd).X = X;
            errorEllipses(structInd).Y = Y;
            errorEllipses(structInd).u = mean(ldaProj);
            structInd = structInd + 1;
        end
    end

%% Get task error ellipses and means 
    taskErrorEllipses = struct('task',[],'posture',[],'X',[],'Y',[],'u',[]);
    structInd = 1;
    for task = taskList
        for posture = postureList
            allTrialAvg = [];
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            for i = 1:size(tempTrajStruct,2)
                allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allZSmoothFR.obs);
            end
            ldaProj = allTrialAvg*taskLDA;
            [X,Y,Z] = error_ellipse_Patrick(ldaProj,pctCI);
            taskErrorEllipses(structInd).task = task;
            taskErrorEllipses(structInd).posture = posture;
            taskErrorEllipses(structInd).X = X;
            taskErrorEllipses(structInd).Y = Y;
            taskErrorEllipses(structInd).u = mean(ldaProj);
            structInd = structInd + 1;
        end
    end


%% Get PTask error ellipses and means for all tasks
    PTErrorEllipses = struct('task',[],'posture',[],'X',[],'Y',[],'Z',[],'u',[]);
    structInd = 1;
    for task = taskList
        for posture = postureList
            allTrialAvg = [];
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            for i = 1:size(tempTrajStruct,2)
                allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allZSmoothFR.obs);
            end
            ldaProj = allTrialAvg*PTaskOrth;
            [X,Y,Z] = error_ellipse_Patrick(ldaProj,pctCI);
            PTErrorEllipses(structInd).task = task;
            PTErrorEllipses(structInd).posture = posture;
            PTErrorEllipses(structInd).X = X;
            PTErrorEllipses(structInd).Y = Y;
            PTErrorEllipses(structInd).Z = Z;
            PTErrorEllipses(structInd).u = mean(ldaProj);
            structInd = structInd + 1;
        end
    end

%% Plot overlapping posture error ellipses
fs = 14;
    figure; hold on;
    for task = taskList
        for posture = postureList
            X = errorEllipses([errorEllipses.task]==task & [errorEllipses.posture]==posture).X; 
            Y = errorEllipses([errorEllipses.task]==task & [errorEllipses.posture]==posture).Y; 
            u = errorEllipses([errorEllipses.task]==task & [errorEllipses.posture]==posture).u; 
            if task == 1
                plot(u(:,1),u(:,2),'.','MarkerSize',27,'Color',pcmap(posture,:));
                plot(X,Y,'Color',pcmap(posture,:),'LineWidth',5);
            elseif task == 2
                plot(u(:,1),u(:,2),'pentagram','MarkerSize',12,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot(X,Y,'--','Color',pcmap(posture,:),'LineWidth',5);
            elseif task == 3
                plot(u(:,1),u(:,2),'square','MarkerSize',12,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot(X,Y,':','Color',pcmap(posture,:),'LineWidth',5);
            end

        end
    end
    ax = gca;
            set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
    ax.TickDir = 'out';
    axis equal;
    xticks([-2,2]); yticks([-2,2])
    postureXLims = ax.XLim; postureYLims = ax.YLim; postureZLims = ax.ZLim;
    if saveFig
       saveas(gcf,fullfile(saveDir,[dataset,'_OverlappingPostureEllipses.svg'])); 
    end

%% Plot scatters of each task individually 
    for task = taskList
        figure; hold on;
        for posture = postureList
            allTrialAvg = [];
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            for i = 1:size(tempTrajStruct,2)
                allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allZSmoothFR.obs);
            end
            ldaProj = allTrialAvg*postureLDA;
            if task ==1 
                plot(ldaProj(:,1),ldaProj(:,2),'.','MarkerSize',20,'Color',pcmap(posture,:))
            elseif task == 2
                plot(ldaProj(:,1),ldaProj(:,2),'pentagram','MarkerSize',10,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            elseif task == 3
                 plot(ldaProj(:,1),ldaProj(:,2),'square','MarkerSize',10,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
            end
        end
        xlim(postureXLims); ylim(postureYLims); zlim(postureZLims);
        ax = gca;
                set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        ax.TickDir = 'out';
        xticks([-2,2]); yticks([-2,2]);
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_task',num2str(task),'_postureScatter.svg'])); 
        end    
    end
    
%% Plot overlapping ellipses, colored by posture, in task space
    figure; hold on;
    for task = taskList
        for posture = postureList
            X = taskErrorEllipses([taskErrorEllipses.task]==task & [taskErrorEllipses.posture]==posture).X; 
            Y = taskErrorEllipses([taskErrorEllipses.task]==task & [taskErrorEllipses.posture]==posture).Y; 
            u = taskErrorEllipses([taskErrorEllipses.task]==task & [taskErrorEllipses.posture]==posture).u; 
            if task == 1
                plot(u(:,1),u(:,2),'.','MarkerSize',25,'Color',pcmap(posture,:));
                plot(X,Y,'Color',pcmap(posture,:),'LineWidth',2);
            elseif task == 2
                plot(u(:,1),u(:,2),'pentagram','MarkerSize',10,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                
                plot(X,Y,'--','Color',pcmap(posture,:),'LineWidth',2);
            elseif task == 3
                plot(u(:,1),u(:,2),'square','MarkerSize',10,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot(X,Y,'-.','Color',pcmap(posture,:),'LineWidth',2);
            end

        end
    end
    ax = gca;
            set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
    ax.TickDir = 'out';
    axis equal;

    if saveFig
       saveas(gcf,fullfile(saveDir,[dataset,'_OverlappingPostureEllipses_TaskSpace.svg'])); 
    end


%% Plot 3d view of covariance ellipses
fs = 14;
    figure; hold on;
    for task = taskList
        for posture = postureList
               X = PTErrorEllipses([PTErrorEllipses.task]==task & [PTErrorEllipses.posture]==posture).X; 
               Y = PTErrorEllipses([PTErrorEllipses.task]==task & [PTErrorEllipses.posture]==posture).Y; 
               Z = PTErrorEllipses([PTErrorEllipses.task]==task & [PTErrorEllipses.posture]==posture).Z; 
               u = PTErrorEllipses([PTErrorEllipses.task]==task & [PTErrorEllipses.posture]==posture).u; 
               colorMat = NaN(size(Z,1),size(Z,2),3);
               for j = 1:size(Z,1)
                   for k = 1:size(Z,2)
                       colorMat(j,k,:) = pcmap(posture,:);
                   end
               end
               h = surf(X,Y,Z,colorMat);
               %set(h, 'FaceAlpha', ellAlpha)
               shading faceted
               %plot3(u(:,1),u(:,2),u(:,3),'.','MarkerSize',20,'Color',pcmap(posture,:));
            
        end       
    end
        xlabel('Posture Dim 1'); ylabel('Posture Dim 2'); zlabel('Task Dim 1')
        ax = gca;
        numTicks = 5;
        ax.XTick = ax.XLim(1):(ax.XLim(2)-ax.XLim(1))/numTicks:ax.XLim(2);
        ax.YTick = ax.YLim(1):(ax.YLim(2)-ax.YLim(1))/numTicks:ax.YLim(2);
        ax.ZTick = ax.ZLim(1):(ax.ZLim(2)-ax.ZLim(1))/numTicks:ax.ZLim(2);
        grid on
        xticklabels({}); yticklabels({}); zticklabels({});
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        ax.ZLim(1) = -7.5
        view([-40 10])
        if saveFig
           saveas(gcf,fullfile(saveDir,[dataset,'_taskXposture3d.svg']));
        end

