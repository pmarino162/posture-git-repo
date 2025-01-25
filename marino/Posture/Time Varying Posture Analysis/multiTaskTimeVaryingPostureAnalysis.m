clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220401 - time-varying posture analysis 1';
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;

%% Load Data
    [Data,~,~,~,~,~,~] = loadEarlData20200314_20211210;

%% Get trajStruct
    binWidth = 25;
    trajFields = {'smoothFR','marker'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
        trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',0},{'state','Step 2','first',0}};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Hold','first',0}};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%% Get posture and target lists by task
    taskList = unique([trajStruct.task]);
    numTasks = size(taskList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);  
    postureList = struct('task',[],'list',[],'num',[]);
    targetList = struct('task',[],'list',[],'num',[]);
    i = 1;
    for task = taskList
        taskTrajStruct = trajStruct([trajStruct.task]==task);
        postureList(i).task = task; 
        postureList(i).list = unique([taskTrajStruct.posture]);
        postureList(i).num = size(postureList(i).list,2);
        targetList(i).task = task; 
        targetList(i).list = unique([taskTrajStruct.target]);
        targetList(i).num = size(targetList(i).list ,2);
        i = i+1;
    end

%% Get minNumTimetamps in condition averages by task
    minNumTimestamps = struct('task',[],'num',[]);
    i = 1;
    for task = taskList
        taskTrajStruct = trajStruct([trajStruct.task]==task);
        numTimestamps = [];
        for j = 1:size(taskTrajStruct,2)
            numTimestamps(j) = length(taskTrajStruct(j).avgSmoothFR.timestamps);
        end
        minNumTimestamps(i).task = task; 
        [minNumTimestamps(i).num,~] = min(numTimestamps);
        i = i+1;
    end

%% Form allTraj containing trial-averaged data for each condition
    allTraj = []; allPos = [];
    for i = 1:size(trajStruct,2)
        task = trajStruct(i).task;
        traj = trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps(task).num,:);
        allTraj = vertcat(allTraj,traj);
        pos = trajStruct(i).avgMarker.traj(1:minNumTimestamps(task).num,1);
        allPos = vertcat(allPos,pos);
    end

%% Do PCA on all data
    [PCA,score,latent,tsquared,explained,mu] = pca(allTraj);     
    
%% Plot x and y hand trajectories for each task
    for task = taskList
        figure
        numTimestamps = minNumTimestamps(task).num;
        for posture = postureList(task).list
            for target = targetList(task).list
                traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj(1:numTimestamps,:); 
                time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps(1:numTimestamps);
                ax1 = subplot(2,1,1)
                    plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                    hold on
                ax2 = subplot(2,1,2)
                    plot(time,traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                    hold on
            end
        end
        subplot(2,1,1)
            ylabel('x (mm)')
        subplot(2,1,2)
            xlabel('time from go (ms)')
            ylabel('y (mm)')
        ax1Range = ax1.YLim(2)-ax1.YLim(1);
        ax2Mid = mean([ax2.YLim(2),ax2.YLim(1)]);
        ax2.YLim = [ax2Mid-ax1Range/2 ax2Mid+ax1Range/2];
        if saveFig
            saveas(gcf,fullfile(saveDir,['Hand Pose Subplot task',num2str(task),'.jpg']))
        end
%             linkaxes([ax1 ax2],'xy')
        figure
        for posture = postureList(task).list
            for target = targetList(task).list
                traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj(1:numTimestamps,:); 
                time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps(1:numTimestamps);
                plot(traj(:,1),traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            end
        end
        xlabel('x (mm)')
        ylabel('y (mm)')
        axis equal
        if saveFig
            saveas(gcf,fullfile(saveDir,['Hand Pose 2d plot task',num2str(task),'.jpg']))
        end
    end
    
%% Do posture LDA on BCI data
    task = 1;
    taskTrajStruct = trajStruct([trajStruct.task]==task);
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = postureList(task).list
       obsStruct(structInd).label = posture;
       allObs = [];
       tempTrajStruct = taskTrajStruct([taskTrajStruct.posture]==posture);
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               traj = tempTrajStruct(i).allSmoothFR(j).traj;
               if size(traj,1) > minNumTimestamps(task).num
                    traj = traj(1:minNumTimestamps(task).num,:);
               end
               allObs = vertcat(allObs,traj);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    [postureLDA] = doLDA(obsStruct);

%% Set up regression struct; collect data for each task into it   
    regStruct = struct('task',[],'X',[],'y',[],'b',[],'proj',[],'projPC',[]);
    structInd = 1;
    for task = 1:4
        %Get taskTrajStruct
        if task == 4
            taskTrajStruct = trajStruct;
        else
            taskTrajStruct = trajStruct([trajStruct.task]==task);
        end
        %Get numBalPts (number of points to keep from each condition to
        %balance across conditions)        
        numBalPts = inf;
        for i = 1:size(taskTrajStruct,2)
            numCondPts = 0;
            numTraj = size(taskTrajStruct(i).allSmoothFR,2);
            for j = 1:numTraj
                traj = taskTrajStruct(i).allSmoothFR(j).traj;
                pos = taskTrajStruct(i).allMarker(j).traj(:,1); 
                if size(pos,1) ~= size(traj,1)
                else
                    numCondPts = numCondPts + size(traj,1);
                end
            end
            if numCondPts < numBalPts
                numBalPts = numCondPts;
            end
        end
        %Get X and y by selecting numBalPts at random from each condition
        X = []; y = [];
        for i = 1:size(taskTrajStruct,2)
            numTraj = size(taskTrajStruct(i).allSmoothFR,2);    
            for j = 1:numTraj
                traj = taskTrajStruct(i).allSmoothFR(j).traj;
                pos = taskTrajStruct(i).allMarker(j).traj(:,1); 
                if size(pos,1) ~= size(traj,1)
                else
                    X = vertcat(X,traj);
                    y = vertcat(y,pos);
                end
            end
            selInd = randsample(size(X,1),numBalPts);
            X = X(selInd,:);
            y = y(selInd,:);
        end
        %Compute b
        kList = [0:5:500];
        resultStruct = struct('k',[],'error',[],'uError',[]);
        n = length(y);
        numFolds = 5;
        c = cvpartition(n,'KFold',numFolds);
        kInd = 1;
        for k = kList
            error = [];
            for fold = 1:numFolds
                idxTrain = training(c,fold);
                idxTest = ~idxTrain;
                b = ridge(y(idxTrain),X(idxTrain,:),k,0);
                yhat = b(1) + X(idxTest,:)*b(2:end); 
                error(fold) = sum((y(idxTest)-yhat).^2);
            end
            resultStruct(kInd).k = k;
            resultStruct(kInd).error = error;
            resultStruct(kInd).uError = mean(error);
            kInd = kInd + 1;
        end
        [~,minErrorInd] = min([resultStruct.uError]);
        k = resultStruct(minErrorInd).k;
        b = ridge(y,X,k,0);
        %Compute proj length
        b = b(2:end)./norm(b(2:end));
        proj = norm(b'*postureLDA);
        projPC = norm(b'*PCA(:,1:2));
        %Fill regStruct
        regStruct(structInd).task = task;
        regStruct(structInd).X = X;
        regStruct(structInd).y = y;
        regStruct(structInd).b = b;
        regStruct(structInd).proj = proj;
        regStruct(structInd).projPC = projPC;
        structInd = structInd + 1;
    end    
    
%% Orthonormalize combinations of axes
    [PPCOrth,~] = qr([postureLDA(:,1),PCA(:,1:2)]); PPCOrth = PPCOrth(:,1:3);
    
%% Plot 3d views of data in BCI posture subspace
    plotTargetList = struct('task',[],'list',[]);
    plotTargetList(1).task = 1; plotTargetList(1).list = [1,5];
    plotTargetList(2).task = 2; plotTargetList(2).list = [1,5];
    plotTargetList(3).task = 3; plotTargetList(3).list = [3,7];
    
    for task = taskList
       figure; hold on;     grid on
       for posture = postureList(task).list
           for target = plotTargetList(task).list
               if any([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target) 
                    traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(1:minNumTimestamps(task).num,:);
                    traj = traj*PPCOrth;
                    plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,1),traj(end,2),traj(end,3),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               end
           end
       end
       xlabel('BCI Posture LDA')
        ylabel('PC 1')
        zlabel('PC 2')
        if saveFig
           saveas(gcf,fullfile(saveDir,['bciPostureAxis3d task',num2str(task),'.fig'])) 
        end
    end

%% Plot 2d views of data in BCI posture subspace
    plotTargetList = struct('task',[],'list',[]);
    plotTargetList(1).task = 1; plotTargetList(1).list = [1,5];
    plotTargetList(2).task = 2; plotTargetList(2).list = [1,5];
    plotTargetList(3).task = 3; plotTargetList(3).list = [3,7];
    
    for task = taskList
       figure; hold on;     grid on
       for posture = postureList(task).list
           for target = plotTargetList(task).list
               if any([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target) 
                    traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(1:minNumTimestamps(task).num,:);
                    traj = traj*postureLDA(:,1:2);
                    plot(traj(:,1),traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                    plot(traj(1,1),traj(1,2),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot(traj(end,1),traj(end,2),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               end
           end
       end
       xlabel('BCI Posture LDA 1')
       ylabel('BCI Posture LDA 2')
        if saveFig
           saveas(gcf,fullfile(saveDir,['bciPostureAxis3d task',num2str(task),'.fig'])) 
        end
    end    
    
%% Plot 2d view 
    plotTargetList = struct('task',[],'list',[]);
    plotTargetList(1).task = 1; plotTargetList(1).list = [1,5];
    plotTargetList(2).task = 2; plotTargetList(2).list = [1,5];
    plotTargetList(3).task = 3; plotTargetList(3).list = [3,7];
    
    for task = taskList
       figure; hold on;     grid on
       for posture = postureList(task).list
           for target = plotTargetList(task).list
               if any([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target) 
                    traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(1:minNumTimestamps(task).num,:);
                    traj = traj*postureLDA(:,1);
                    time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(1:minNumTimestamps(task).num);
                    subplot(2,1,1)
                        plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                        hold on
                    traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj(1:minNumTimestamps(task).num,:);
                    time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps(1:minNumTimestamps(task).num);
                    subplot(2,1,2)
                        plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                        hold on
               end
           end
       end
       %Compute p
       X = regStruct(task).X;
       y = regStruct(task).y;
       yhat = X*postureLDA(:,1);
       p = corrcoef(y,yhat);
       subplot(2,1,1)
            ax = gca;
            xlim = ax.XLim;
            ylim = ax.YLim;
            title(['p = ',num2str(p(1,2))])
            %text(xlim(1),ylim(2)-,['p = ',num2str(p(1,2))])
            ylabel('BCI Posture LDA')
       subplot(2,1,2)
            xlabel('time (ms)')
            ylabel('x (mm)')   
       if saveFig
           saveas(gcf,fullfile(saveDir,['bciPostureAxis2d task',num2str(task),'.jpg'])) 
       end
    end

    
%% Compute pMat to see how well one posture axis predicts other postures
    pMat = zeros(4,4);
    for trainTask = 1:4
        b = regStruct(trainTask).b;
        for testTask = 1:4
            xTest = regStruct(testTask).X;
            yTest = regStruct(testTask).y;
            yhat = xTest*b; 
            p = corrcoef(yTest,yhat);
            pMat(trainTask,testTask) = p(1,2);
        end
    end

    figure
    imagesc(pMat)
    c = colorbar;
    caxis([0 1])
    c.Label.String = '\rho';
    
    xticks([1:4])
    xticklabels({'BCI','Reaching','Iso','All'})
    xlabel('test task')
    
    yticks([1:4])
    yticklabels({'BCI','Reaching','Iso','All'})
    ylabel('train task')
    
    if saveFig
        saveas(gcf,fullfile(saveDir,'pMatrix.jpg')) 
    end
        
%% Plot 2d view for cross-projections
    plotTargetList = struct('task',[],'list',[]);
    plotTargetList(1).task = 1; plotTargetList(1).list = [1,5];
    plotTargetList(2).task = 2; plotTargetList(2).list = [1,5];
    plotTargetList(3).task = 3; plotTargetList(3).list = [3,7];
    
    for trainTask = 1:4
        for testTask = taskList
           figure; hold on;     grid on
           for posture = postureList(testTask).list
               for target = plotTargetList(testTask).list
                   if any([trajStruct.task]==testTask & [trajStruct.posture]==posture & [trajStruct.target]==target) 
                        traj = trajStruct([trajStruct.task]==testTask & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(1:minNumTimestamps(testTask).num,:);
                        b = regStruct(trainTask).b;
                        traj = traj*b;
                        time = trajStruct([trajStruct.task]==testTask & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(1:minNumTimestamps(testTask).num);
                        subplot(2,1,1)
                            plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                            hold on
                        traj = trajStruct([trajStruct.task]==testTask & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj(1:minNumTimestamps(testTask).num,:);
                        time = trajStruct([trajStruct.task]==testTask & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps(1:minNumTimestamps(testTask).num);
                        subplot(2,1,2)
                            plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                            hold on
                   end
               end
           end
           %Compute p
           X = regStruct(testTask).X;
           y = regStruct(testTask).y;
           yhat = X*b;
           p = corrcoef(y,yhat);
           subplot(2,1,1)
                ax = gca;
                xlim = ax.XLim;
                ylim = ax.YLim;
                title(['p = ',num2str(p(1,2))])
                %text(xlim(1),ylim(2)-,['p = ',num2str(p(1,2))])
                ylabel('Posture regression')
           subplot(2,1,2)
                xlabel('time (ms)')
                ylabel('x (mm)')   
           if saveFig
                saveas(gcf,fullfile(saveDir,['postureReg train',num2str(trainTask),' test',num2str(testTask),'.jpg'])) 
           end
        end
    end


%% Plot projection lengths
    figure
    proj = [regStruct.proj];
    projPC = [regStruct.projPC];
    bar([proj;projPC]')
    ylabel('projection length')
    ax = gca;
    ax.YLim = [0 1];
    legend('Proj Into Posture LDA','Proj Into PC')
    xticklabels({'BCI','Reaching','Iso','All'})
    xlabel('train task')
    
%% Visualize performance on averages
    trainTask = 2;
    testTask = 2;
    
    figure; hold on
    for posture = [1,3,5]
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==testTask & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj; 
            b = regStruct(trainTask).b;
            traj = traj*postureLDA(:,1);
            %traj = traj*b;
            time = trajStruct([trajStruct.task]==testTask & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
            subplot(3,1,1)
                plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            traj = trajStruct([trajStruct.task]==testTask & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj; 
            time = trajStruct([trajStruct.task]==testTask & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps;
            subplot(3,1,2)
                plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            subplot(3,1,3)
                plot(time,traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
        end
    end
    subplot(3,1,2)
        ylabel('x (mm)')
    subplot(3,1,3)
        xlabel('time (ms)')
        ylabel('y (mm)')       


%% Do ridge regression to predict position from neural activity
    task = 2;
    taskTrajStruct = trajStruct([trajStruct.task]==task);
    X = []; y = [];
    for i = 1:size(taskTrajStruct,2)
        numTraj = size(taskTrajStruct(i).allSmoothFR,2);
        for j = 1:numTraj
            numTimestamps = length(taskTrajStruct(i).allSmoothFR(j).timestamps);
            traj = taskTrajStruct(i).allSmoothFR(j).traj;
            pos = taskTrajStruct(i).allMarker(j).traj(:,1); 
            if size(pos,1) ~= size(traj,1)
            else
                numTimestamps = size(pos,1);
                if numTimestamps > minNumTimestamps(task).num
                    traj = traj(1:minNumTimestamps(task).num,:);
                    pos = pos(1:minNumTimestamps(task).num,:);
                end
                X = vertcat(X,traj);
                y = vertcat(y,pos);
            end
        end
    end

    kList = [0,1,5,10:5:50,75,100,200];
    resultStruct = struct('k',[],'error',[],'uError',[]);
    
    n = length(y);
    numFolds = 5;
    c = cvpartition(n,'KFold',numFolds);
    kInd = 1;
    for k = kList
        error = [];
        for fold = 1:numFolds
            idxTrain = training(c,fold);
            idxTest = ~idxTrain;
            b = ridge(y(idxTrain),X(idxTrain,:),k,0);
            yhat = b(1) + X(idxTest,:)*b(2:end); 
            error(fold) = sum((y(idxTest)-yhat).^2);
        end
        resultStruct(kInd).k = k;
        resultStruct(kInd).error = error;
        resultStruct(kInd).uError = mean(error);
        kInd = kInd + 1;
    end

    figure
    hold on
    for i = 1:size(resultStruct,2)
        k = resultStruct(i).k;
        error = resultStruct(i).error;
        uError = resultStruct(i).uError;
        semilogx(k,uError,'.','MarkerSize',10)
    end
    semilogx([resultStruct.k],[resultStruct.uError])
    xlabel('k')
    ylabel('error')
    
    figure
    k = 40;
    b = ridge(y,X,k,0);
    yhat = b(1) + X(idxTest,:)*b(2:end);  
    scatter(y(idxTest),yhat)
    hold on
    plot(y(idxTest),y(idxTest))
    xlabel('Actual x Pos')
    ylabel('Predicted x Pos')
    hold off
    postureReg = b(2:end);
    
%% Do ridge regression on subsets of data - see how each predicts others
    pMat = zeros(4,4);
    k = 40;
    pStruct = struct('trainTask',[],'testTask',[],'p',[],'up',[]);
    structInd = 1;
    for trainTask = 1:4
        xTrain = []; yTrain = [];
        if trainTask == 4
            trainTaskTrajStruct = trajStruct;
        else
            trainTaskTrajStruct = trajStruct([trajStruct.task]==trainTask);
        end
        for i = 1:size(trainTaskTrajStruct,2)
            numTraj = size(trainTaskTrajStruct(i).allSmoothFR,2);
            tempTask = trainTaskTrajStruct(i).task;
            for j = 1:numTraj
                numTimestamps = length(trainTaskTrajStruct(i).allSmoothFR(j).timestamps);
                traj = trainTaskTrajStruct(i).allSmoothFR(j).traj;
                pos = trainTaskTrajStruct(i).allMarker(j).traj(:,1); 
                if size(pos,1) ~= size(traj,1)
                else
                    if numTimestamps > minNumTimestamps(tempTask).num
                        traj = traj(1:minNumTimestamps(tempTask).num,:);
                        pos = pos(1:minNumTimestamps(tempTask).num,:);
                    end
                    xTrain = vertcat(xTrain,traj);
                    yTrain = vertcat(yTrain,pos);
                end
            end
        end
        b = ridge(yTrain,xTrain,k,0);
        for testTask = 1:4
            testTask
            xTest = []; yTest = [];
            if testTask == 4
                testTaskTrajStruct = trajStruct;
            else
                testTaskTrajStruct = trajStruct([trajStruct.task]==testTask);
            end
            for i = 1:size(testTaskTrajStruct,2)
                numTraj = size(testTaskTrajStruct(i).allSmoothFR,2);
                tempTask = testTaskTrajStruct(i).task;
                for j = 1:numTraj
                    numTimestamps = length(testTaskTrajStruct(i).allSmoothFR(j).timestamps);
                    traj = testTaskTrajStruct(i).allSmoothFR(j).traj;
                    pos = testTaskTrajStruct(i).allMarker(j).traj(:,1); 
                    if size(pos,1) ~= size(traj,1)
                    else
                        numTimestamps = size(pos,1);
                        if numTimestamps > minNumTimestamps(tempTask).num
                            traj = traj(1:minNumTimestamps(tempTask).num,:);
                            pos = pos(1:minNumTimestamps(tempTask).num,:);
                        end
                        xTest = vertcat(xTest,traj);
                        yTest = vertcat(yTest,pos);
                    end
                end
            end
            %pStruct(structInd).trainTask = trainTask;
            %pStruct(structInd).testTask = testTask;
            %pStruct(structInd).p = p;
            %pStruct(structInd).up = up;
            %structInd = structInd + 1;
            yhat = b(1) + xTest*b(2:end); 
            p = corrcoef(yTest,yhat);
            pMat(trainTask,testTask) = p(1,2);
        end
    end
    
    figure
    imagesc(pMat)
    
    
    

    
%% Add projections to trajStruct
    for i = 1:size(trajStruct,2)
        %trajStruct(i).PCA.traj = trajStruct(i).avgSmoothFR.traj*PCA;
        %trajStruct(i).PPCOrth.traj = trajStruct(i).avgSmoothFR.traj*PPCOrth;
        %trajStruct(i).PRPCOrth.traj = trajStruct(i).avgSmoothFR.traj*PRPCOrth;
        trajStruct(i).postureLDA.traj = trajStruct(i).avgSmoothFR.traj*postureLDA;
        trajStruct(i).postureReg.traj = trajStruct(i).avgSmoothFR.traj*postureReg;
    end
    
%% Look at reaching position trajectories to understand time
    figure; hold on
    task = 2;
    for posture = [1,3,5]
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).postureReg.traj; 
            time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
            subplot(3,1,1)
                plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj; 
            time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps;
            subplot(3,1,2)
                plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            subplot(3,1,3)
                plot(time,traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
        end
    end
    subplot(3,1,2)
        ylabel('x (mm)')
    subplot(3,1,3)
        xlabel('time (ms)')
        ylabel('y (mm)')       
    
    
    
    
%% Plot marker pose vs time
    figure; hold on
    task = 2;
    timePts = 1:18;
    plotTimePts = timePts;
    numTimestamps = 18;
    for posture = postureList(task).list
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(timePts,:); 
            traj = b(1) + traj*b(2:end);  
            plot(plotTimePts,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
            plotTimePts = plotTimePts + numTimestamps;
        end
    end
    
%     figure; hold on
    timeInd = 1;
    task = 2;
    timePts = 1:18;
    plotTimePts = timePts;
    numTimestamps = 18;
    for posture = postureList
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj(timePts,:); 
%             time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(timePts,:);
            plot(plotTimePts,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
            plotTimePts = plotTimePts + numTimestamps;
        end
    end
       
    
    
%% Look at reaching position trajectories to understand time
    figure; hold on
    task = 2;
    for posture = 1
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj; 
            time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps;
            subplot(2,1,1)
                plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            subplot(2,1,2)
                plot(time,traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
        end
    end
    subplot(2,1,1)
        ylabel('x (mm)')
    subplot(2,1,2)
        xlabel('time (ms)')
        ylabel('y (mm)')
    




%% Look at reaching position trajectories to understand time
    figure; hold on
    task = 2;
    for posture = [1,3,5]
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).postureReg.traj; 
            time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
            subplot(3,1,1)
                plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj; 
            time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps;
            subplot(3,1,2)
                plot(time,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
            subplot(3,1,3)
                plot(time,traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                hold on
        end
    end
    subplot(3,1,2)
        ylabel('x (mm)')
    subplot(3,1,3)
        xlabel('time (ms)')
        ylabel('y (mm)')    

%% Plot marker pose vs time
    figure; hold on
    timeInd = 1;
    task = 2;
    timePts = 1:18;
    plotTimePts = timePts;
    numTimestamps = 18;
    for posture = postureList
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).PPCOrth.traj(timePts,:); 
%             time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(timePts,:);
            plot(plotTimePts,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
            plotTimePts = plotTimePts + numTimestamps;
        end
    end
    
%     figure; hold on
    timeInd = 1;
    task = 2;
    timePts = 1:18;
    plotTimePts = timePts;
    numTimestamps = 18;
    for posture = postureList
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj(timePts,:); 
%             time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(timePts,:);
            plot(plotTimePts,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
            plotTimePts = plotTimePts + numTimestamps;
        end
    end
    
    %% Plot marker pose vs time
    figure; hold on
    timeInd = 1;
    task = 2;
    timePts = 1:18;
    plotTimePts = timePts;
    numTimestamps = 18;
    for posture = postureList
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(timePts,:); 
            traj = b(1) + traj*b(2:end);  
            plot(plotTimePts,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
            plotTimePts = plotTimePts + numTimestamps;
        end
    end
    
%     figure; hold on
    timeInd = 1;
    task = 2;
    timePts = 1:18;
    plotTimePts = timePts;
    numTimestamps = 18;
    for posture = postureList
        for target = [1,3,5]
            traj = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj(timePts,:); 
%             time = trajStruct([trajStruct.task]==task & [trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(timePts,:);
            plot(plotTimePts,traj(:,1),'Color',pcmap(posture,:),'LineWidth',2);
            plotTimePts = plotTimePts + numTimestamps;
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