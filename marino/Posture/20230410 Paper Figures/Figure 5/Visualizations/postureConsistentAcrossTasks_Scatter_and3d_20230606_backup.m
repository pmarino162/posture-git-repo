clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\Consistent Posture Signal';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    pcmap = pcmap([1,3,5],:);

%% Parameters
   numManPCs = 10;     %Num PCs to project data into before analysis
    pctCI = 95;
    ellAlpha = 0.4;
    
%% Load Data, subselect
    dataset = 'E20200314';
    [Data,zScoreParams] = loadData(dataset);
    
    %Get trajStruct
    binWidth = 25; kernelStdDev = 25;
    trialInclStates = struct('trialName','','inclStates',[]);
    trajFields = {'zSmoothFR','markerVel','marker'};
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    

    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
        trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',0},{'state','Step 2','first',0}};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-300},{'kin','moveOnsetTime','first',-100}};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);


    taskIDs = struct('ID',[],'task','');
    taskIDs(1).ID = 1; taskIDs(1).task = 'BC';
    taskIDs(2).ID = 2; taskIDs(2).task = 'HC';
    taskIDs(3).ID = 3; taskIDs(3).task = 'Iso';
    taskIDs(4).ID = 4; taskIDs(4).task = 'All';
    
%% For each timecourse in trajstruct, get single point
    for i = 1:size(trajStruct,2)
        %Individual trials
        for j = 1:size(trajStruct(i).allSmoothFR,2)
           trajStruct(i).allSmoothFR(j).trialAvg = mean(trajStruct(i).allSmoothFR(j).traj); 
        end
    end

%% Reduce dimensionality of observations
        allObs = [];
        for i = 1:size(trajStruct,2)
            for j = 1:size(trajStruct(i).allSmoothFR,2)
                obs = trajStruct(i).allSmoothFR(j).trialAvg;
                allObs = vertcat(allObs,obs);
            end
        end
        [coeff,score,latent,tsquared,explained,mu] = pca(allObs); 
        for i = 1:size(trajStruct,2)
            for j = 1:size(trajStruct(i).allSmoothFR,2)
                trajStruct(i).allSmoothFR(j).trialAvg = trajStruct(i).allSmoothFR(j).trialAvg*coeff(:,1:numManPCs);
            end
        end
        
%% Get posture, target, and task lists
    postureList = unique([trajStruct.posture]);
    targetList = unique([trajStruct.target]);
    taskList = unique([trajStruct.task]);
    
 %% Create taskObsStruct, balancing observations by target within each task x posture
    %Get minimum number of trials to any target direction within each
    %task x posture
    minNumObsStruct = struct('task',[],'posture',[],'minNumObs',[]);
    structInd = 1;
    for task = 1:3
        for posture = 1:3
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            tempMinNumCondObs = 999999999;
            for i = 1:size(tempTrajStruct,2)
                numCondObs = size(vertcat(tempTrajStruct(i).allSmoothFR.trialAvg),1);
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

    %Build taskObsStruct
    taskObsStruct = struct('task',[],'obs',[],'labels',[]);
    structInd = 1;
    for task = 1:4
        obs = []; labels = [];
        if task < 4
            for posture = 1:3
                tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
                numObsToUsePerTarget = minNumObsStruct([minNumObsStruct.task]==task & [minNumObsStruct.posture]==posture).minNumObs;
                for i = 1:size(tempTrajStruct,2)
                    condObs = vertcat(tempTrajStruct(i).allSmoothFR.trialAvg);
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


%% Do LDA by posture and task on all trial avgerages
    obs = taskObsStruct(4).obs; labels = taskObsStruct(4).labels;
    numClasses = 3;
    LDAproj = fisherLDA(obs, labels);
    [LDAproj,~] = qr(LDAproj);
    LDAproj = LDAproj(:,1:numClasses-1);
    postureLDA = LDAproj;
%     %Flip axes so that P1 is on the right and P2 is on the top 
%     P1Proj = mean(obsStruct([obsStruct.label]==1).allObs)*postureLDA;
%     P2Proj = mean(obsStruct([obsStruct.label]==2).allObs)*postureLDA;
%     if P1Proj(:,1) < P2Proj(:,1)
%         postureLDA(:,1) = -1.*postureLDA(:,1);
%     end
%     if P1Proj(:,2) > P2Proj(:,2)
%         postureLDA(:,2) = -1.*postureLDA(:,2);
%     end
%   
    
%% Do LDA by posture on all trial averages 
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = postureList
       obsStruct(structInd).label = posture;
       allObs = [];
       tempTrajStruct = trajStruct([trajStruct.posture]==posture);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               trialAvg = tempTrajStruct(i).allSmoothFR(j).trialAvg;
               allObs = vertcat(allObs,trialAvg);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    [postureLDA] = doLDA(obsStruct);
    
    %Flip axes so that P1 is on the right and P2 is on the top 
    P1Proj = mean(obsStruct([obsStruct.label]==1).allObs)*postureLDA;
    P2Proj = mean(obsStruct([obsStruct.label]==2).allObs)*postureLDA;
    if P1Proj(:,1) < P2Proj(:,1)
        postureLDA(:,1) = -1.*postureLDA(:,1);
    end
    if P1Proj(:,2) > P2Proj(:,2)
        postureLDA(:,2) = -1.*postureLDA(:,2);
    end
    
%% Do task LDA on trial averages; create PTaskOrth
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for task = 1:3
       obsStruct(structInd).label = task;
       allObs = [];
       tempTrajStruct = trajStruct([trajStruct.task]==task);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               trialAvg = tempTrajStruct(i).allSmoothFR(j).trialAvg;
               allObs = vertcat(allObs,trialAvg);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    [taskLDA] = doLDA(obsStruct);
    [PTaskOrth,~] = qr([postureLDA,taskLDA(:,1)]); PTaskOrth = PTaskOrth(:,1:3);
    %Flip axes so that P1 is on the right and P3 is on the top 
    P1Proj = mean(obsStruct([obsStruct.label]==1).allObs)*PTaskOrth;
    P3Proj = mean(obsStruct([obsStruct.label]==3).allObs)*PTaskOrth;
    if P1Proj(:,1) > P3Proj(:,1)
        PTaskOrth(:,1) = -1.*PTaskOrth(:,1);
    end
    if P1Proj(:,2) > P3Proj(:,2)
        PTaskOrth(:,2) = -1.*PTaskOrth(:,2);
    end


%% Get posture error ellipses and means for all tasks
    errorEllipses = struct('task',[],'posture',[],'X',[],'Y',[],'u',[]);
    structInd = 1;
    for task = taskList
        for posture = postureList
            allTrialAvg = [];
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            for i = 1:size(tempTrajStruct,2)
                allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allSmoothFR.trialAvg);
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

%% Get PTask error ellipses and means for all tasks
    PTErrorEllipses = struct('task',[],'posture',[],'X',[],'Y',[],'Z',[],'u',[]);
    structInd = 1;
    for task = taskList
        for posture = postureList
            allTrialAvg = [];
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            for i = 1:size(tempTrajStruct,2)
                allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allSmoothFR.trialAvg);
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
    figure; hold on;
    for task = taskList
        for posture = postureList
            X = errorEllipses([errorEllipses.task]==task & [errorEllipses.posture]==posture).X; 
            Y = errorEllipses([errorEllipses.task]==task & [errorEllipses.posture]==posture).Y; 
            u = errorEllipses([errorEllipses.task]==task & [errorEllipses.posture]==posture).u; 
            if task == 1
                plot(u(:,1),u(:,2),'.','MarkerSize',25,'Color',pcmap(posture,:));
                plot(X,Y,'Color',pcmap(posture,:),'LineWidth',2);
            elseif task == 2
                plot(u(:,1),u(:,2),'square','MarkerSize',10,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot(X,Y,'--','Color',pcmap(posture,:),'LineWidth',2);
            elseif task == 3
                plot(u(:,1),u(:,2),'pentagram','MarkerSize',10,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                plot(X,Y,'-.','Color',pcmap(posture,:),'LineWidth',2);
            end

        end
    end
    ax = gca;
    ax.TickDir = 'out';
    ax.XTick = [];
    ax.YTick = [];
    
%% Plot scatters of each task individually 
    for task = taskList
        figure; hold on;
        for posture = postureList
            allTrialAvg = [];
            tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
            for i = 1:size(tempTrajStruct,2)
                allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allSmoothFR.trialAvg);
            end
            ldaProj = allTrialAvg*postureLDA;
            plot(ldaProj(:,1),ldaProj(:,2),'.','MarkerSize',20,'Color',pcmap(posture,:))
        end
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
               set(h, 'FaceAlpha', ellAlpha)
               shading interp
               plot3(u(:,1),u(:,2),u(:,3),'.','MarkerSize',20,'Color',pcmap(posture,:));
            
        end       
    end
            xlabel('Posture Dim 1'); ylabel('Posture Dim 2'); zlabel('Task Dim 1')
        xticklabels({}); yticklabels({}); zticklabels({});
        grid on
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
       % xlim(postureXLimits); ylim(postureYLimits);
        view([0 90])
        if saveFig
           saveas(gcf,fullfile(saveDir,[dataset,'_taskXposture2d.svg']));
           saveas(gcf,fullfile(saveDir,[dataset,'_taskXposture2d.fig']));
        end
        view([-40 10])
        if saveFig
           saveas(gcf,fullfile(saveDir,[dataset,'_taskXposture3d.svg']));
           saveas(gcf,fullfile(saveDir,[dataset,'_taskXposture3d.fig']));
        end


% %% Local function for performing LDA
%      function [LDAproj] = doLDA(obsStruct)
%         %Preallocate
%         minNumObs = min([obsStruct.numObs]);
%         numClasses = size(obsStruct,2);
%         numDims = size(obsStruct(1).allObs,2);
%         obs = NaN(minNumObs*numClasses,numDims); labels =  NaN(minNumObs*numClasses,1); 
%         %Fill
%         k = 1;
%         for i = 1:size(obsStruct,2)
%             totalClassObs = size(obsStruct(i).allObs,1);
%             obs(k:k+minNumObs-1,:) = obsStruct(i).allObs(randsample(totalClassObs,minNumObs),:);
%             labels(k:k+minNumObs-1,1) = ones(minNumObs,1).*obsStruct(i).label;
%             k = k+minNumObs;
%         end
%         LDAproj = fisherLDA(obs, labels);
%         [LDAproj,~] = qr(LDAproj);
%         LDAproj = LDAproj(:,1:numClasses-1);
%      end  
     
%% Local function for performing LDA - no balancing 
     function [LDAproj] = doLDA(obsStruct)
        numClasses = size(obsStruct,2);
        numDims = size(obsStruct(1).allObs,2);
        obs = []; labels = [];
        for i = 1:numel(obsStruct)
            obs = vertcat(obs,obsStruct(i).allObs);
            labels = vertcat(labels,obsStruct(i).label*ones(size(obsStruct(i).allObs,1),1));
        end
        LDAproj = fisherLDA(obs, labels);
        [LDAproj,~] = qr(LDAproj);
        LDAproj = LDAproj(:,1:numClasses-1);
     end 