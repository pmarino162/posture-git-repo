clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\Consistent Posture Signal';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%% Set pctCI and ellAlpha for error ellipses
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
        trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',25},{'state','Step 2','first',0}};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-300},{'kin','moveOnsetTime','first',-50}};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {{'state','Target','first',75},{'state','Target','first',200}};
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

    taskIDs = struct('ID',[],'task','');
    taskIDs(1).ID = 1; taskIDs(1).task = 'BC';
    taskIDs(2).ID = 2; taskIDs(2).task = 'HC';
    taskIDs(3).ID = 3; taskIDs(3).task = 'Iso';
    
%% For each timecourse in trajstruct, get single point
    for i = 1:size(trajStruct,2)
        %Individual trials
        for j = 1:size(trajStruct(i).allSmoothFR,2)
           trajStruct(i).allSmoothFR(j).trialAvg = mean(trajStruct(i).allSmoothFR(j).traj); 
        end
    end

%% Get posture, target, and task lists
    postureList = unique([trajStruct.posture]);
    targetList = unique([trajStruct.target]);
    taskList = unique([trajStruct.task]);
    
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
    
    %Flip axes so that P1 is on the right and P3 is on the top 
    P1Proj = mean(obsStruct([obsStruct.label]==1).allObs)*postureLDA;
    P3Proj = mean(obsStruct([obsStruct.label]==3).allObs)*postureLDA;
    if P1Proj(:,1) < P3Proj(:,1)
        postureLDA(:,1) = -1.*postureLDA(:,1);
    end
    if P1Proj(:,2) > P3Proj(:,2)
        postureLDA(:,2) = -1.*postureLDA(:,2);
    end
    
%% For each task, do targetLDA
    targetLDAStruct = struct('task',[],'targetLDA',[]);
    ldaStructInd = 1;
    for task = 1:3
        obsStruct = struct('label',[],'numObs',[],'allObs',[]);
        obsStructInd = 1;
        taskTrajStruct = trajStruct([trajStruct.task]==task);
        taskTargetList = unique([taskTrajStruct.target]);
        for target = taskTargetList
           allObs = [];
           tempTrajStruct = taskTrajStruct([taskTrajStruct.target]==target & [taskTrajStruct.task]==task);
           for i = 1:size(tempTrajStruct,2)
               for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
                   trialAvg = tempTrajStruct(i).allSmoothFR(j).trialAvg;
                   allObs = vertcat(allObs,trialAvg);
               end
           end
           numObs = size(allObs,1);
           obsStruct(obsStructInd).numObs = numObs;
           obsStruct(obsStructInd).allObs = allObs;
           obsStruct(obsStructInd).label = target;
           obsStructInd = obsStructInd + 1; 
        end
        [targetLDA] = doLDA(obsStruct);
        
        %Flip axes so that T1 is on the right and T3 is on the top
        if task == 1 || task == 2
            T1Proj = mean(obsStruct([obsStruct.label]==1).allObs)*targetLDA;
            T3Proj = mean(obsStruct([obsStruct.label]==3).allObs)*targetLDA;
            if T1Proj(:,1) < T3Proj(:,1)
                targetLDA(:,1) = -1.*targetLDA(:,1);
            end
            if T1Proj(:,2) > T3Proj(:,2)
                targetLDA(:,2) = -1.*targetLDA(:,2);
            end
        elseif task == 3
            T3Proj = mean(obsStruct([obsStruct.label]==3).allObs)*targetLDA;
            T7Proj = mean(obsStruct([obsStruct.label]==7).allObs)*targetLDA;
            if T7Proj(:,1) > T3Proj(:,1)
                targetLDA(:,1) = -1.*targetLDA(:,1);
            end
        end
        targetLDAStruct(ldaStructInd).task = task;
        targetLDAStruct(ldaStructInd).targetLDA = targetLDA;
        ldaStructInd = ldaStructInd + 1;   
    end
    

%% Get posture error ellipses for all tasks
    condGroup = struct('taskList',[],'postureList',[],'targetList',[]);
    condGroup(1).taskList = [1]; condGroup(1).postureList = [1,3,5]; condGroup(1).targetList = 1:8;
    condGroup(2).taskList = [2]; condGroup(2).postureList = [1,3,5]; condGroup(2).targetList = 1:8;
    condGroup(3).taskList = [3]; condGroup(3).postureList = [1,3,5]; condGroup(3).targetList = 1:8;
    
    errorEllipses = struct('task',[],'posture',[],'X',[],'Y',[]);
    structInd = 1;
    figure
    for condGroupInd = 1:3
        taskList = condGroup(condGroupInd).taskList;
        postureList = condGroup(condGroupInd).postureList;
        for task = taskList
            for posture = postureList
                if any([trajStruct.task]==task & [trajStruct.posture]==posture)
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
                    structInd = structInd + 1;
                    fill(X,Y,pcmap(posture,:),'FaceAlpha',ellAlpha);
                    hold on
                end
            end
        end 
    end
    
%% Plot all posture points and error ellipses to get posture space axis limits
    fs = 14;
    condGroup = struct('taskList',[],'postureList',[],'targetList',[]);
    condGroup(1).taskList = [1:3]; condGroup(1).postureList = [1,3,5]; condGroup(1).targetList = 1:8;
    for condGroupInd = 1
        taskList = condGroup(condGroupInd).taskList;
        postureList = condGroup(condGroupInd).postureList;
        targetList = condGroup(condGroupInd).targetList;
        figure;
        hold on
        for task = taskList
            for posture = postureList
                X = errorEllipses([errorEllipses.task]==task & [errorEllipses.posture]==posture).X;
                Y = errorEllipses([errorEllipses.task]==task & [errorEllipses.posture]==posture).Y;
                fill(X,Y,pcmap(posture,:),'FaceAlpha',ellAlpha);
                for target = targetList
                    if any([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target)
                        allTrialAvg = [];
                        tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target);
                        for i = 1:size(tempTrajStruct,2)
                            allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allSmoothFR.trialAvg);
                        end
                        ldaProj = allTrialAvg*postureLDA;
                        plot(ldaProj(:,1),ldaProj(:,2),'.','MarkerSize',20,'Color',pcmap(posture,:));
                    end
                end
            end
        end
        xticklabels({}); yticklabels({}); zticklabels({}); 
        grid on
        xlabel('Posture Dim 1')
        ylabel('Posture Dim 2')
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        ax = gca; 
        postureXLimits = ax.XLim;
        postureYLimits = ax.YLim;
        postureZLimits = ax.ZLim;
    end

%% Plot each task in posture LDA space separately 
    fs = 14;
    
    condGroup = struct('taskList',[],'postureList',[],'targetList',[]);
    condGroup(1).taskList = [1:3]; condGroup(1).postureList = [1,3,5]; condGroup(1).targetList = 1:8;
    condGroup(2).taskList = [1]; condGroup(2).postureList = [1]; condGroup(2).targetList = 3;
    condGroup(3).taskList = [1]; condGroup(3).postureList = [1]; condGroup(3).targetList = 3;
    condGroup(4).taskList = [1]; condGroup(4).postureList = [1]; condGroup(4).targetList = 1:8; 
    condGroup(5).taskList = [1]; condGroup(5).postureList = [1,3]; condGroup(5).targetList = 1:8;
    condGroup(6).taskList = [1]; condGroup(6).postureList = [1,3,5]; condGroup(6).targetList = 1:8;
    condGroup(7).taskList = [1]; condGroup(7).postureList = [1,3,5]; condGroup(7).targetList = 1:8;
    condGroup(8).taskList = [1]; condGroup(8).postureList = []; condGroup(8).targetList = [];
    condGroup(9).taskList = [3]; condGroup(9).postureList = [1]; condGroup(9).targetList = 3;
    condGroup(10).taskList = [3]; condGroup(10).postureList = [1]; condGroup(10).targetList = [3,7];
    condGroup(11).taskList = [3]; condGroup(11).postureList = [1,3,5]; condGroup(11).targetList = [3,7];
    condGroup(12).taskList = [1,3]; condGroup(12).postureList = []; condGroup(12).targetList = [];
    condGroup(13).taskList = [1:3]; condGroup(13).postureList = [1,3,5]; condGroup(13).targetList = [1:8];
    condGroup(14).taskList = [1:3]; condGroup(14).postureList = []; condGroup(14).targetList = [];
    
    for condGroupInd = 1:14
        taskList = condGroup(condGroupInd).taskList;
        postureList = condGroup(condGroupInd).postureList;
        targetList = condGroup(condGroupInd).targetList;
        figure;
        hold on
        if condGroupInd >=7 %Starting w CG7, include BCI posture CI ellipse
            for posture = [1,3,5]
                X = errorEllipses([errorEllipses.task]==1 & [errorEllipses.posture]==posture).X;
                Y = errorEllipses([errorEllipses.task]==1 & [errorEllipses.posture]==posture).Y;
                fill(X,Y,pcmap(posture,:),'FaceAlpha',ellAlpha);
                if condGroupInd >=12 %Starting w CG12, also include Iso posture CI ellipse
                    X = errorEllipses([errorEllipses.task]==3 & [errorEllipses.posture]==posture).X;
                    Y = errorEllipses([errorEllipses.task]==3 & [errorEllipses.posture]==posture).Y;
                    fill(X,Y,pcmap(posture,:),'FaceAlpha',ellAlpha);
                end
                if condGroupInd >=14 %Starting w CG12, also include Reach posture CI ellipse
                    X = errorEllipses([errorEllipses.task]==2 & [errorEllipses.posture]==posture).X;
                    Y = errorEllipses([errorEllipses.task]==2 & [errorEllipses.posture]==posture).Y;
                    fill(X,Y,pcmap(posture,:),'FaceAlpha',ellAlpha);
                end    
            end
        end
        for task = taskList
            for posture = postureList
                for target = targetList
                    if any([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target)
                        allTrialAvg = [];
                        tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target);
                        for i = 1:size(tempTrajStruct,2)
                            allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allSmoothFR.trialAvg);
                        end
                        ldaProj = allTrialAvg*postureLDA;
                        if condGroupInd == 2
                            plot(ldaProj(1,1),ldaProj(1,2),'.','MarkerSize',20,'Color',pcmap(posture,:));
                        else 
                            plot(ldaProj(:,1),ldaProj(:,2),'.','MarkerSize',20,'Color',pcmap(posture,:));
                        end
                    end
                end
            end
        end
        
        xlim(postureXLimits)
        ylim(postureYLimits)
        zlim(postureZLimits)
        xticklabels({}); yticklabels({}); zticklabels({}); 
        grid on
        xlabel('Posture Dim 1')
        ylabel('Posture Dim 2')
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        if saveFig           
           if condGroupInd == 7
              saveas(gcf,fullfile(saveDir,[dataset,'_task',num2str(taskList),'_posture',num2str(postureList),'_target',num2str(targetList),'_2dPostureLDA_wEllipse.svg']));
           elseif condGroupInd == 2
              saveas(gcf,fullfile(saveDir,[dataset,'_task',num2str(taskList),'_posture',num2str(postureList),'_target',num2str(targetList),'_2dPostureLDA_OnePoint.svg']));
           else
              saveas(gcf,fullfile(saveDir,[dataset,'_task',num2str(taskList),'_posture',num2str(postureList),'_target',num2str(targetList),'_2dPostureLDA.svg']));
           end
        end
    end

%% Plot each task in target LDA space separately. First plot all data for each task, and get x and y limits. 
    %Then plot subsets of data from each task in that space.
    fs = 14;
    taskAxisLims = struct('task',[],'xlimits',[],'ylimits',[]);
    isoScatter = struct('posture',[],'target',[],'scatter',[]);
    isoScatterStructInd = 1;
    
    condGroup = struct('taskList',[],'postureList',[],'targetList',[]);
    condGroup(1).taskList = [1]; condGroup(1).postureList = [1,3,5]; condGroup(1).targetList = 1:8;
    condGroup(2).taskList = [2]; condGroup(2).postureList = [1,3,5]; condGroup(2).targetList = 1:8;
    condGroup(3).taskList = [3]; condGroup(3).postureList = [1,3,5]; condGroup(3).targetList = 1:8;
    
    condGroup(4).taskList = [1]; condGroup(4).postureList = [1]; condGroup(4).targetList = 3;
    condGroup(5).taskList = [1]; condGroup(5).postureList = [1]; condGroup(5).targetList = 3;
    condGroup(6).taskList = [1]; condGroup(6).postureList = [1]; condGroup(6).targetList = 1:8;
    condGroup(7).taskList = [1]; condGroup(7).postureList = [1,3]; condGroup(7).targetList = 1:8;
    condGroup(8).taskList = [3]; condGroup(8).postureList = [1]; condGroup(8).targetList = 3;
    condGroup(9).taskList = [3]; condGroup(9).postureList = [1]; condGroup(9).targetList = [3,7];
    condGroup(10).taskList = [3]; condGroup(10).postureList = [1,3,5]; condGroup(10).targetList = [3,7];
    condGroup(11).taskList = [2]; condGroup(11).postureList = [1,3,5]; condGroup(11).targetList = [1:8];

    for condGroupInd =1:11
        taskList = condGroup(condGroupInd).taskList;
        postureList = condGroup(condGroupInd).postureList;
        targetList = condGroup(condGroupInd).targetList;
        figure;
        hold on
        for task = taskList
            targetLDA = targetLDAStruct([targetLDAStruct.task]==task).targetLDA;
            for posture = postureList
                for target = targetList
                    if any([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target)
                        allTrialAvg = [];
                        tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target);
                        for i = 1:size(tempTrajStruct,2)
                            allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allSmoothFR.trialAvg);
                        end
                        ldaProj = allTrialAvg*targetLDA;
                        if condGroupInd == 3
                            isoScatter(isoScatterStructInd).posture = posture;
                            isoScatter(isoScatterStructInd).target = target;
                            isoScatter(isoScatterStructInd).scatter = 0.15*randn(size(ldaProj,1),1);
                            isoScatterStructInd = isoScatterStructInd + 1;
                        end
                        if task == 3
                            scatterPts = isoScatter([isoScatter.posture]==posture & [isoScatter.target]==target).scatter;
                            plot(scatterPts,ldaProj(:,1),'.','MarkerSize',20,'Color',tcmap(target,:));
                        elseif condGroupInd == 4
                            plot(ldaProj(1,1),ldaProj(1,2),'.','MarkerSize',20,'Color',tcmap(target,:));
                        else
                            plot(ldaProj(:,1),ldaProj(:,2),'.','MarkerSize',20,'Color',tcmap(target,:));
                        end
                    end
                end
            end
        end
        
        if ismember(condGroupInd,[1:3])
            ax = gca; 
            task = taskList;
            taskAxisLims(task).task = task;
            taskAxisLims(task).xlimits = ax.XLim;
            taskAxisLims(task).ylimits = ax.YLim;
            taskAxisLims(task).zlimits = ax.ZLim;  
        end
        xlim(taskAxisLims(task).xlimits)
        ylim(taskAxisLims(task).ylimits)
        zlim(taskAxisLims(task).zlimits)
        xlabel('Target Dim 1')
        ylabel('Target Dim 2')
        if task == 3
            xlabel('')
            ylabel('Target Dim 1')
            xlim(taskAxisLims(task).ylimits)
        end
        xticklabels({}); yticklabels({}); zticklabels({}); 
        grid on
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
       
        if saveFig
            if condGroupInd == 4
                saveas(gcf,fullfile(saveDir,[dataset,'_task',num2str(taskList),'_posture',num2str(postureList),'_target',num2str(targetList),'_2dTargetLDA_OnePoint.svg']));
            else
                saveas(gcf,fullfile(saveDir,[dataset,'_task',num2str(taskList),'_posture',num2str(postureList),'_target',num2str(targetList),'_2dTargetLDA.svg']));
            end
        end
    end
  
   
%% Empty posture space
    figure
    p =plot([],'.','MarkerSize',20,'Color',[0 0 0]);
    xlim(postureXLimits); ylim(postureYLimits); zlim(postureZLimits)
    xticklabels({}); yticklabels({}); zticklabels({}); 
    grid on
    xlabel('Posture Dim 1'); ylabel('Posture Dim 2')
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    if saveFig
       saveas(gcf,fullfile(saveDir,[dataset,'EmptyPostureSpace.svg']));
    end

%% Empty BCI target space
    figure
    task = 1;
    p =plot([],'.','MarkerSize',20,'Color',[0 0 0]);
    xlim(taskAxisLims([taskAxisLims.task]==task).xlimits); ylim(taskAxisLims([taskAxisLims.task]==task).ylimits); zlim(taskAxisLims([taskAxisLims.task]==task).zlimits)
    xticklabels({}); yticklabels({}); zticklabels({}); 
    grid on
    xlabel('Target Dim 1'); ylabel('Target Dim 2')
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    if saveFig
       saveas(gcf,fullfile(saveDir,[dataset,'EmptyBCITargetSpace.svg']));
    end    
    
%% Empty Iso target space
    figure
    task = 3;
    p =plot([],'.','MarkerSize',20,'Color',[0 0 0]);
    xlim(taskAxisLims([taskAxisLims.task]==task).xlimits); ylim(taskAxisLims([taskAxisLims.task]==task).ylimits); zlim(taskAxisLims([taskAxisLims.task]==task).zlimits)
    xticklabels({}); yticklabels({}); zticklabels({}); 
    grid on
    xlabel(''); ylabel('Target Dim 1')
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    if saveFig
       saveas(gcf,fullfile(saveDir,[dataset,'EmptyIsoTargetSpace.svg']));
    end   
    
%% Do LDA by task on all trial averages; get PTaskOrth and make sure postures are oriented correctly
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
    
%% Plot task x posture means
    condGroup = struct('taskList',[],'postureList',[]);
    condGroup(1).taskList = [1:3]; condGroup(1).postureList = [1,3,5]; 
    %condGroup(2).taskList = [1:3]; condGroup(2).postureList = [1,3,5]; 
    
    for condGroupInd = 1
        taskList = condGroup(condGroupInd).taskList;
        postureList = condGroup(condGroupInd).postureList;
        figure
        hold on
        for task = taskList
            for posture = postureList
                allTrialAvg = [];
                tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
                for i = 1:size(tempTrajStruct,2)
                    allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allSmoothFR.trialAvg);
                end
                ldaProj = allTrialAvg*PTaskOrth;
                mu = mean(ldaProj);
                [X,Y,Z] = error_ellipse_Patrick(ldaProj,pctCI);
                colorMat = NaN(size(Z,1),size(Z,2),3);
               for j = 1:size(Z,1)
                   for k = 1:size(Z,2)
                       colorMat(j,k,:) = pcmap(posture,:);
                   end
               end
               h = surf(X,Y,Z,colorMat);
               set(h, 'FaceAlpha', ellAlpha)
               shading interp
               plot3(mu(:,1),mu(:,2),mu(:,3),'.','MarkerSize',20,'Color',pcmap(posture,:));
%            plot3(ldaProj(:,1),ldaProj(:,2),ldaProj(:,3),'.','MarkerSize',20,'Color',pcmap(posture,:));
            end
        end
        xlabel('Posture Dim 1'); ylabel('Posture Dim 2'); zlabel('Task Dim 1')
        xticklabels({}); yticklabels({}); zticklabels({});
        grid on
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        xlim(postureXLimits); ylim(postureYLimits);
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
    end
    ax = gca;
    xlimits = ax.XLim; ylimits = ax.YLim; zlimits = ax.ZLim;



%% Plot all tasks in 3d space; rotate upward to show task axis
    fileName = fullfile(saveDir,[dataset,'_RotatingUpToShowTask.mp4']);
    f = figure; set(gcf,'color','white'); hold on
    %daspect([1 1 1])
%     f.OuterPosition = [0 0 1000 1000];
    %f.Position = [10 10 500 1000];
    %ax = gca
    %ax.ActivePositionProperty = 'position'
    for task = [1:3]
        for posture = [1,3,5]
           allTrialAvg = [];
           tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
           for i = 1:size(tempTrajStruct,2)
                allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allSmoothFR.trialAvg);
           end
           ldaProj = allTrialAvg*PTaskOrth;
           mu = mean(ldaProj);
           [X,Y,Z] = error_ellipse_Patrick(ldaProj,pctCI);
           colorMat = NaN(size(Z,1),size(Z,2),3);
           for j = 1:size(Z,1)
               for k = 1:size(Z,2)
                   colorMat(j,k,:) = pcmap(posture,:);
               end
           end
           h = surf(X,Y,Z,colorMat);
           set(h, 'FaceAlpha', ellAlpha)
           shading interp
           plot3(mu(:,1),mu(:,2),mu(:,3),'.','MarkerSize',20,'Color',pcmap(posture,:));
%            plot3(ldaProj(:,1),ldaProj(:,2),ldaProj(:,3),'.','MarkerSize',20,'Color',pcmap(posture,:));
        end
    end
    %xlabel('Posture Dim 1'); ylabel('Posture Dim 2'); zlabel('Task Dim 1')
    xticklabels({}); yticklabels({}); zticklabels({});
    
    grid on
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
    %axis vis3d
    ax = gca; ax.XLim = postureXLimits; ax.YLim = postureXLimits; ax.ZLim = zlimits;
        OptionZ.FrameRate=12;OptionZ.Duration=5.5;OptionZ.Periodic=false;
    viewVec = [0,90; -40,10];
    CaptureFigVid(viewVec,fileName,OptionZ)


% %% Fade from BCI to empty space
%     f = figure; set(gcf,'color','white'); hold on
%     fileName = fullfile(saveDir,[dataset,'_BCItoBlankFade.mp4']);
%     fps = 5;
%     v = VideoWriter(fileName); v.FrameRate = fps; open(v);
%     
%     condGroupInd = 5;
%     taskList = condGroup(condGroupInd).taskList; postureList = condGroup(condGroupInd).postureList; targetList = condGroup(condGroupInd).targetList;
%     for frame = 1:10
%         for task = taskList
%             for posture = postureList
%                     if any([trajStruct.task]==task & [trajStruct.posture]==posture)
%                         allTrialAvg = [];
%                         tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
%                         for i = 1:size(tempTrajStruct,2)
%                             allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allSmoothFR.trialAvg);
%                         end
%                         ldaProj = allTrialAvg*postureLDA;
%                     end
%                     color = pcmap(posture,:);
%                     color = interp1([0 10],[pcmap(posture,:); 1 1 1],frame);
%         
%                     p =plot(ldaProj(:,1),ldaProj(:,2),'.','MarkerSize',20,'Color',color);
%             end
%         end
%         
%         if frame == 10
%             clf
%             p =plot([],'.','MarkerSize',20,'Color',color);
%                         xlim(xlimits); ylim(ylimits); zlim(zlimits)
%             xticklabels({}); yticklabels({}); zticklabels({}); 
%             grid on
%             xlabel('Posture Dim 1'); ylabel('Posture Dim 2')
%             set(gca,'fontname','arial'); set(gca,'fontsize',fs)
%         end
%         
%         if frame == 1
%             xlim(xlimits); ylim(ylimits); zlim(zlimits)
%             xticklabels({}); yticklabels({}); zticklabels({}); 
%             grid on
%             xlabel('Posture Dim 1'); ylabel('Posture Dim 2')
%             set(gca,'fontname','arial'); set(gca,'fontsize',fs)
%         end
%         M(frame) = getframe(gcf); writeVideo(v,M(frame));
%     end
%     close(v);
%     movie(M,1,fps)
% %% Plot iso force ellipse, then add in BCI
%     [PTaskOrth,~] = qr([postureLDA,taskLDA(:,1)]); PTaskOrth = PTaskOrth(:,1:3);
%     
%     condGroup = struct('taskList',[],'postureList',[],'targetList',[]);
%     condGroup(1).taskList = [3]; condGroup(1).postureList = [1,3,5]; condGroup(1).targetList = 1:8;
%     condGroup(2).taskList = [1,3]; condGroup(2).postureList = [1,3,5]; condGroup(2).targetList = 1:8;
%     
%         
%     for condGroupInd = 1:2
%         taskList = condGroup(condGroupInd).taskList;
%         postureList = condGroup(condGroupInd).postureList;
%         figure
%         hold on
%         for task = taskList
%             for posture = postureList
%                 allTrialAvg = [];
%                 tempTrajStruct = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture);
%                 for i = 1:size(tempTrajStruct,2)
%                     allTrialAvg = vertcat(allTrialAvg,tempTrajStruct(i).allSmoothFR.trialAvg);
%                 end
%                 ldaProj = allTrialAvg*PTaskOrth;
%                 mu = mean(ldaProj);
%                 %ell = error_ellipse(ldaProj);
%                 %ell2 = error_ellipse(ldaProj(:,2:3));
%                 xr = 25; yr = 25; zr = 25;
%                 [X,Y,Z] = ellipsoid(mu(:,1),mu(:,2),mu(:,3),xr,yr,zr);
%                 colorMat = NaN(size(Z,1),size(Z,2),3);
%                for j = 1:size(Z,1)
%                    for k = 1:size(Z,2)
%                        colorMat(j,k,:) = pcmap(posture,:);
%                    end
%                end
%                h = surf(X,Y,Z,colorMat);
%                set(h, 'FaceAlpha', 0.2)
%                shading interp
%                plot3(mu(:,1),mu(:,2),mu(:,3),'.','MarkerSize',20,'Color',pcmap(posture,:));
% %            plot3(ldaProj(:,1),ldaProj(:,2),ldaProj(:,3),'.','MarkerSize',20,'Color',pcmap(posture,:));
%             end
%         end
%         xlabel('Posture Dim 1'); ylabel('Posture Dim 2'); zlabel('Task Dim 1')
%         xticklabels({}); yticklabels({}); zticklabels({});
%         grid on
%         set(gca,'fontname','arial')
%         set(gca,'fontsize',fs)
%         ax = gca;
%         ax.XLim = xlimits; ax.YLim = ylimits; ax.ZLim = zlimits;
%         view([0 90])
%         if saveFig
%            saveas(gcf,fullfile(saveDir,[dataset,'condGroup',num2str(condGroupInd),'_taskXposture2d.svg']));
%            %saveas(gcf,fullfile(saveDir,[dataset,'_taskXposture2d.fig']));
%         end
% 
%     end
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