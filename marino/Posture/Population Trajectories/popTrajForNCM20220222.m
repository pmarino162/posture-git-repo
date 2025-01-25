clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
%     task = 'iso';
%     task  = 'BCI';
%     task = 'planning';
    task = 'reaching';
%     task = 'multijoint BCI';
    switch task
        case 'BCI'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\BCI Figures';
        case 'planning'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\Planning Figures';
        case 'iso'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\Iso Figures';
        case 'reaching'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\Reaching Figures';
        case 'multijoint BCI'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\Multijoint BCI Figures';
    end
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;
    
%% Load Data, subselect, get traj struct
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;

    switch task
        case 'BCI'
        [Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
        [Data] = subselectForTrajDist(Data,task);
        trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
        condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
        %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-150},{'state','Step 2','first',0}};
        trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
        trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

        case 'multijoint BCI'
        [Data] = loadEarlData20210901_20211210;
        taskID = [Data.conditionData];
        taskID = [taskID.taskID];
        BCData = Data(taskID==1);
        HCData = Data(taskID==2);
        [BCData] = subselectForTrajDist(BCData,'BCI');
        Data = BCData; clearvars BCData;
        %[HCData] = subselectForTrajDist(HCData,'reaching');
        trialInclStates(1).trialName = {'BCI Center Out'};
        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
        trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
   
        case 'iso'
        [Data] = loadEarlData20200116_20211210();
        [Data] = subselectForTrajDist(Data,task);
        trialInclStates(1).trialName = {'IsometricForce_1D'};
        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%         trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
        trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
        trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    
        case 'reaching'
        [Data] = loadEarlData20210706_20211210; 
        [Data] = subselectForTrajDist(Data,task);
        trialInclStates(1).trialName = {'GridReaching'};
        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
        trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
        fourPostureTrajStruct = trajStruct(ismember([trajStruct.posture],[1:4]));
    %     trajStruct = trajStruct(ismember([trajStruct.posture],[1,2,3,5,6,7]));

        case 'planning'
        [Data] = loadEarlData20210706_20211210; 
        [Data] = subselectForMarginalization(Data,task);
         trialInclStates(1).trialName = {'GridReaching'};
        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        trialInclStates(1).inclStates = {{'state','Delay','first',0},{'state','Delay','first',250}};
        trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
        trajStruct = trajStruct(ismember([trajStruct.posture],[1:4]));
    end
    
%% Get timestamps dist and min number timestamps
    numTimestamps = [];
    numCondTraj = [];
    for i = 1:size(trajStruct,2)
       numTraj = size(trajStruct(i).allSmoothFR,2);
       numCondTraj = [numCondTraj,numTraj];
       for j = 1:numTraj
          numTimestamps = [numTimestamps,size(trajStruct(i).allSmoothFR(j).timestamps,2)]; 
       end
    end
    figure
    histogram(numTimestamps)
    xlabel('Number of 25ms bins')
    ylabel('Number of trials')
    
    figure
    histogram(numCondTraj)
    xlabel('Number of trials')
    ylabel('Number of conditions')
    
% Get minimum number of trials and timestamps
    [minNumTimestamps,i] = min(numTimestamps);
    [minNumCondTraj,i] = min(numCondTraj);
    
%% Get posture and target lists
    trajStruct = trajStruct;
    switch task
        case {'BCI','planning','iso','multijoint BCI'}
            postureList = unique([trajStruct.posture]);
        case 'reaching'
            postureList = unique([fourPostureTrajStruct.posture]);
    end
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
      
%% Do PCA on condition averages for visualization; get CIs and Var expl
    allAvgs = NaN(minNumTimestamps*size(trajStruct,2),size(trajStruct(1).allSmoothFR(1).traj,2));
    j = 1;
    for i = 1:size(trajStruct,2)
        allAvgs(j:j+minNumTimestamps-1,:) = trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps,:);
        j = j+minNumTimestamps;
    end
    [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);       
    
%% Peform LDA by posture on all data 
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = postureList
       obsStruct(structInd).label = posture;
       allObs = [];
       switch task
           case {'BCI','multijoint BCI'}
               tempTrajStruct = trajStruct([trajStruct.posture]==posture);
           case 'planning'
               tempTrajStruct = trajStruct([trajStruct.posture]==posture);
           case 'reaching'
               tempTrajStruct = fourPostureTrajStruct([fourPostureTrajStruct.posture]==posture);
           case 'iso'
               tempTrajStruct = trajStruct([trajStruct.posture]==posture);
       end
       
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
               traj = tempTrajStruct(i).allSmoothFR(j).traj;
               traj = traj(1:minNumTimestamps,:);
               allObs = vertcat(allObs,traj);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    
    [postureLDA] = doLDA(obsStruct);

  
%% Peform LDA by target on all data - around 250ms after cue onset
    %Get number of points for each condition
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;
    for target = targetList
       obsStruct(structInd).label = target;
       allObs = [];
       switch task
           case {'BCI','multijoint BCI'}
                tempTrajStruct = trajStruct([trajStruct.target]==target);
                timeToUse = 125;
                numOnEitherSide = 1;
           case 'iso'
                tempTrajStruct = trajStruct([trajStruct.target]==target);
                timeToUse = 175;
                numOnEitherSide = 1;
           case 'reaching'
                tempTrajStruct = trajStruct([trajStruct.target]==target);
                timeToUse = 100;
                numOnEitherSide = 2;
           case 'planning'
                tempTrajStruct = trajStruct([trajStruct.target]==target);
                timeToUse = 175;
                numOnEitherSide = 1;
       end
       for i = 1:size(tempTrajStruct,2)
           for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
               timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
               traj = tempTrajStruct(i).allSmoothFR(j).traj;
               [~,timeToUseInd] = min(abs(timestamps-timeToUse));
               traj = traj(timeToUseInd-numOnEitherSide:timeToUseInd+numOnEitherSide,:);
               allObs = vertcat(allObs,traj);
           end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    
    [targetLDA] = doLDA(obsStruct);
    

%% Get Post Targ Orth
    switch task
        case {'BCI','reaching','planning','multijoint BCI'}
        [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
        postTargOrth = postTargOrth(:,1:3);
        case 'iso'
        [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1)]);  
        postTargOrth = postTargOrth(:,1:2);
    end
    
    
%% Do PCA in null space of posture space 
    postureAxisNullBasis = null(postureLDA(:,1)');
    nullProj = allAvgs*postureAxisNullBasis;
    [postureAxisNullBasisPC,score,latent,tsquared,explained,nullMu] = pca(nullProj);
    postureAxisNullBasis = postureAxisNullBasis*postureAxisNullBasisPC;
        
    pxtNullBasis = null(postTargOrth');
    nullProj = allAvgs*pxtNullBasis;
    [pxtNullBasisPC,score,latent,tsquared,explained,nullMu] = pca(nullProj);
    pxtNullBasis = pxtNullBasis*pxtNullBasisPC;
    
%% Get angles
    acosd(dot(postureLDA(:,1),targetLDA(:,1)))
    acosd(dot(postureLDA(:,1),targetLDA(:,2)))
    
%% Add Projections to trajStruct; do some stats 
    totalVar = trace(cov(allAvgs));
    for i = 1:size(trajStruct,2)
        %Add average traces to trajStruct
        trajStruct(i).avgPCA.traj = (trajStruct(i).avgSmoothFR.traj-mu)*coeff;
        trajStruct(i).avgPostureLDA.traj = trajStruct(i).avgSmoothFR.traj*postureLDA;
        trajStruct(i).avgTargetLDA.traj = trajStruct(i).avgSmoothFR.traj*targetLDA;
        trajStruct(i).avgPCApNull.traj = trajStruct(i).avgSmoothFR.traj*postureAxisNullBasis;
        trajStruct(i).avgPTOrth.traj = trajStruct(i).avgSmoothFR.traj*postTargOrth;
        trajStruct(i).avgPCApxtNull.traj = trajStruct(i).avgSmoothFR.traj*pxtNullBasis;
        %Get VAF
        trajStruct(i).avgPCA.VAF = 100.*(diag(cov(allAvgs*coeff))')./totalVar;
        trajStruct(i).avgPostureLDA.VAF = 100.*(diag(cov(allAvgs*postureLDA))')./totalVar;
        trajStruct(i).avgTargetLDA.VAF =  100.*(diag(cov(allAvgs*targetLDA))')./totalVar;
        trajStruct(i).avgPCApNull.VAF =  100.*(diag(cov(allAvgs*postureAxisNullBasis))')./totalVar;
        trajStruct(i).avgPTOrth.VAF =  100.*(diag(cov(allAvgs*postTargOrth))')./totalVar;
        trajStruct(i).avgPCApxtNull.VAF = 100.*(diag(cov(allAvgs*pxtNullBasis))')./totalVar;
        %Add all traces to trajStruct; Compute 95% CI
        numTrials = size(trajStruct(i).allSmoothFR,2);
        allPCA = NaN(minNumTimestamps,size(trajStruct(i).avgPCA.traj,2),numTrials);
        allPostureLDA = NaN(minNumTimestamps,size(trajStruct(i).avgPostureLDA.traj,2),numTrials);
        allTargetLDA = NaN(minNumTimestamps,size(trajStruct(i).avgTargetLDA.traj,2),numTrials);
        allPCApNull = NaN(minNumTimestamps,size(trajStruct(i).avgPCApNull.traj,2),numTrials);
        allPTOrth = NaN(minNumTimestamps,size(trajStruct(i).avgPTOrth.traj,2),numTrials);
        allPCApxtNull = NaN(minNumTimestamps,size(trajStruct(i).avgPCApxtNull.traj,2),numTrials);
        for j = 1:numTrials
          allPCA(:,:,j) = (trajStruct(i).allSmoothFR(j).traj(1:minNumTimestamps,:)-mu)*coeff;
          allPostureLDA(:,:,j) = trajStruct(i).allSmoothFR(j).traj(1:minNumTimestamps,:)*postureLDA;
          allTargetLDA(:,:,j) = trajStruct(i).allSmoothFR(j).traj(1:minNumTimestamps,:)*targetLDA;
          allPCApNull(:,:,j) = trajStruct(i).allSmoothFR(j).traj(1:minNumTimestamps,:)*postureAxisNullBasis;
          allPTOrth(:,:,j) = trajStruct(i).allSmoothFR(j).traj(1:minNumTimestamps,:)*postTargOrth;
          allPCApxtNull(:,:,j) = trajStruct(i).allSmoothFR(j).traj(1:minNumTimestamps,:)*pxtNullBasis;
        end
        %Get confidence intervals
        trajStruct(i).avgPCA.CI = 1.96.*std(allPCA,0,3)/sqrt(numTrials);
        trajStruct(i).avgPostureLDA.CI = 1.96.*std(allPostureLDA,0,3)/sqrt(numTrials);
        trajStruct(i).avgTargetLDA.CI = 1.96.*std(allTargetLDA,0,3)/sqrt(numTrials);
        trajStruct(i).avgPCApNull.CI = 1.96.*std(allPCApNull,0,3)/sqrt(numTrials);
        trajStruct(i).avgPTOrth.CI = 1.96.*std(allPTOrth,0,3)/sqrt(numTrials);
        trajStruct(i).avgPCApxtNull.CI = 1.96.*std(allPCApxtNull,0,3)/sqrt(numTrials);
    end
    
%% Top 3 PC's
    switch task
        case 'BCI'
            plotPostureList = [1,5];
            plotTargetList = [1,3,5,7];
            timePts = 1:8;
        case 'multijoint BCI'
            plotPostureList = [4,5];
            plotTargetList = [1,3,5];
            timePts = 1:12;
        case 'iso'
%             plotPostureList = [1,5];
%             plotTargetList = [3,7];
            plotPostureList = [1:5];
            plotTargetList = [3,7];
            timePts = 1:10;
        case 'reaching'
            plotPostureList = [2,7];
            plotTargetList = [5,3];
            timePts = 1:15;
%             plotPostureList = [1,2,3,4];
%             plotTargetList = [4,5,6];
%             timePts = 1:15;

        case 'planning'
            plotPostureList = [1,4];
            plotTargetList = [1,5];
            timePts = 1:8;
    end
%     pcmap = customSummer;
%     pcmap(6,:) = customSummer(5,:);
%     pcmap(7,:) = customSummer(5,:);
    figure
    hold on
     for posture = plotPostureList
       for target = plotTargetList
           if any([trajStruct.posture]==posture & [trajStruct.target]==target)
               traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPCA.traj(timePts,:); 
    %            traj = traj-[0,traj(1,2:3)];
               plot3(traj(:,1),traj(:,2),traj(:,4),'Color',pcmap(posture,:),'LineWidth',2);
               plot3(traj(1,1),traj(1,2),traj(1,4),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
               plot3(traj(end,1),traj(end,2),traj(end,4),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
           end
       end
     end
     
    VAF = round(trajStruct(1).avgPCA.VAF);
    xlabel(['PC 1 (',num2str(VAF(1)),'%)'])
    ylabel(['PC 2 (',num2str(VAF(2)),'%)'])
    zlabel(['PC 4 (',num2str(VAF(4)),'%)']) 
    grid on
    axis equal
        
    if saveFig
       saveas(gcf,fullfile(saveDir,'Top3PC.fig')); 
    end
    
%     set(gcf, 'Color', [1 1 1])
%     vidFileName = fullfile(saveDir,'Top3PCVid');
%     OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
%     viewVec = [90,10; 0,10; -90,10; -180,10; -270,10];
%     CaptureFigVid(viewVec, vidFileName,OptionZ)


    
%% Posture LDA and top 2 PCs
    switch task
        case 'BCI'
            plotPostureList = 1:5;
            plotTargetList = [1,3,5,7];
            timePts = 1:8;
        case 'multijoint BCI'
            plotPostureList = [4,5];
            plotTargetList = [1,3,5];
            timePts = 1:12;
        case 'iso'
            plotPostureList = 1:5;
            plotTargetList = [3,7];
            timePts = 1:8;
        case 'reaching'
            plotPostureList = [1,4];
            plotTargetList = [2,6];
            timePts = 1:10;
        case 'planning'
            plotPostureList = [1,4];
            plotTargetList = [1,5];
            timePts = 1:8;
    end
    
    figure
    hold on
     for posture = plotPostureList
       for target = plotTargetList
           postureLDA = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPostureLDA.traj; 
           PCApNull = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPCApNull.traj;
           traj = [postureLDA(timePts,1),PCApNull(timePts,[1,2])];
%            traj = traj-[0,traj(1,2:3)];
           plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
           plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
           plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
       end
    end
    postureVAF = trajStruct(1).avgPostureLDA.VAF;
    PCApNullVAF = trajStruct(1).avgPCApNull.VAF;
    xlabel(['Posture LDA 1','(',num2str(round(postureVAF(1))),'%)'])
    ylabel(['PC 1 ','(',num2str(round(PCApNullVAF(1))),'%)'])
    zlabel(['PC 2 ','(',num2str(round(PCApNullVAF(2))),'%)']) 
    axis equal
    grid on
    
%% PxT Projection 
pcmap = customSummer;
    switch task
        case 'BCI'
            plotPostureList = 1:5;
            plotTargetList = [1,3,5,7];
            timePts = 1:8;
        case 'multijoint BCI'
            plotPostureList = [4,5];
            plotTargetList = [1,3,5];
            timePts = 1:12;
        case 'iso'
            plotPostureList = 1:5;
            plotTargetList = [3,7];
            timePts = 1:8;
        case 'reaching'
            plotPostureList = [1:4];
            plotTargetList = [1:8];
            timePts = 1:10;
        case 'planning'
            plotPostureList = [1:4];
            plotTargetList = [1:8];
            timePts = 1:8;
    end

    figure
    hold on
    switch task
        case {'BCI','reaching','planning','multijoint BCI'}
             for posture = plotPostureList
               for target = plotTargetList
                   if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                   traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPTOrth.traj(timePts,:); 
                   plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',1.5);
                   plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',25,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                   plot3(traj(end,1),traj(end,2),traj(end,3),'o','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',[1 1 1]);
                   end
               end
             end
            pxtVAF = trajStruct(1).avgPTOrth.VAF;
            xlabel(['Posture LDA 1','(',num2str(round(pxtVAF(1))),'%)'])
            ylabel(['Target LDA 1 ','(',num2str(round(pxtVAF(2))),'%)'])
            zlabel(['Target LDA 2 ','(',num2str(round(pxtVAF(3))),'%)']) 
        case 'iso'
            for posture = plotPostureList
               for target = plotTargetList
                   traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPTOrth.traj(timePts,:); 
                   plot(traj(:,1),traj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
                   plot(traj(1,1),traj(1,2),'.','MarkerSize',25,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                   plot(traj(end,1),traj(end,2),'o','MarkerSize',10,'Color',pcmap(posture,:),'MarkerFaceColor',[1 1 1]);
               end
            end
            pxtVAF = trajStruct(1).avgPTOrth.VAF;
            xlabel(['Posture LDA 1','(',num2str(round(pxtVAF(1))),'%)'])
            ylabel(['Target LDA 1 ','(',num2str(round(pxtVAF(2))),'%)'])
    end
%     axis equal
    grid on
    
    if saveFig
       saveas(gcf,fullfile(saveDir,'PxT.fig')); 
    end
    
    
%     set(gcf, 'Color', [1 1 1])
%     vidFileName = fullfile(saveDir,'PxTVid');
%     OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
%     viewVec = [90,10; 0,10; -90,10; -180,10; -270,10];
%     CaptureFigVid(viewVec, vidFileName,OptionZ)
    
%% PxT and Null Timecourses
    switch task
        case 'BCI'
            plotPostureList = [1,5];
            plotTargetList = [3,7];
            timePts = 1:8;
    end
    
    figure
    for target = plotTargetList
        for posture = plotPostureList
            if ~isempty(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target))
               time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(timePts); 
               pxtTraj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPTOrth.traj(timePts,:); 
               pxtNullTraj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPCApxtNull.traj(timePts,:); 
               traj = [pxtTraj(:,1:3),pxtNullTraj];
               pxtCI = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPTOrth.CI(timePts,:);
               pxtNullCI = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPCApxtNull.CI(timePts,:);
               CI = [pxtCI(:,1:3),pxtNullCI(:,1:7)];
               for i = 1:10
                   subplot(2,5,i)
                   shadedErrorBar(time,traj(:,i),CI(:,i),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                   hold on
               end
            end
        end
%     pxtVAF = 
%     pxtNullVAF = 
%     VAF = 
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
           ylabel('Posture Dim 1')
       elseif plotInd ==2
           ylabel('Target Dim 1')
       elseif plotInd == 3
           ylabel('Target Dim 2')
       else
           ylabel(['Null ',num2str(plotInd-3)])
       end
    end

    end

     

%% PCA Timecourses
    switch task
        case 'BCI'
            plotPostureList = [1,5];
            plotTargetList = [1,5];
            timePts = 1:8;
        case 'multijoint BCI'
            plotPostureList = [4,5];
            plotTargetList = [1,3,5];
            timePts = 1:12;
        case 'iso'
            plotPostureList = [1,5];
            plotTargetList = [3,7];
            timePts = 1:8;
         case 'reaching'
            plotPostureList = [2,3];
            plotTargetList = [7];
            timePts = 1:15;
         case 'planning'
            plotPostureList = [1,4];
            plotTargetList = [1,5];
            timePts = 1:8;
    end
    
    figure
    for target = plotTargetList
        for posture = plotPostureList
            if ~isempty(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target))
               time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(timePts); 
               traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPCA.traj(timePts,:); 
               CI = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPCA.CI(timePts,:);
               for i = 1:10
                   subplot(2,5,i)
                   shadedErrorBar(time,traj(:,i),CI(:,i),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)})
                   hold on
               end
            end
        end
        
    VAF = round(trajStruct(1).avgPCA.VAF);
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
       ylabel(['PC ',num2str(plotInd),' (',num2str(VAF(plotInd)),'%)'])
    end

    end

%      
%     if saveFig
%        saveas(gcf,fullfile(saveDir,task,['P',num2str(cond1P),'T',num2str(cond1T),'vs','P',num2str(cond2P),'T',num2str(cond2T),'TrajComparisonShifted.jpg'])); 
%     end
%     close all

%% Visualize some averages
    switch task
        case 'BCI'
            timePts = 1:8;
            cond1P = 1;
            cond1T = 1;
            cond2P = 5;
            cond2T = 1;
        case 'multijoint BCI'
            timePts = 1:12;
            cond1P = 4;
            cond1T = 7;
            cond2P = 5;
            cond2T = 7;
        case 'iso'
            timePts = 1:11;
            cond1P = 1;
            cond1T = 3;
            cond2P = 5;
            cond2T = 3;
        case 'reaching'
            timePts = 1:15;
            cond1P = 2;
            cond1T = 5;
            cond2P = 7;
            cond2T = 5;
        case 'planning'
            plotPostureList = [1,4];
            timePts = 1:8;
            cond1P = 1;
            cond1T = 5;
            cond2P = 4;
            cond2T = 5;
    end




    cond1Traj = trajStruct([trajStruct.posture]==cond1P & [trajStruct.target]==cond1T).avgPCA.traj(timePts,:);
    cond1Time = trajStruct([trajStruct.posture]==cond1P & [trajStruct.target]==cond1T).avgSmoothFR.timestamps(timePts);
    cond1CI =  trajStruct([trajStruct.posture]==cond1P & [trajStruct.target]==cond1T).avgPCA.CI(timePts,:);
    
    cond2Traj = trajStruct([trajStruct.posture]==cond2P & [trajStruct.target]==cond2T).avgPCA.traj(timePts,:);
    cond2Time = trajStruct([trajStruct.posture]==cond2P & [trajStruct.target]==cond2T).avgSmoothFR.timestamps(timePts);
    cond2CI = trajStruct([trajStruct.posture]==cond2P & [trajStruct.target]==cond2T).avgPCA.CI(timePts,:);
    
    alpha = getOptimalAlpha(cond1Traj,cond2Traj,length(timePts));
    
    cond2ShiftTraj = cond2Traj + alpha;
    
    %Unshifted
    f = figure; f.Position = [100 100 1500 500];
    for i = 1:10
        subplot(2,5,i)
        shadedErrorBar(cond1Time,cond1Traj(:,i),cond1CI(:,i),'lineprops',{'LineWidth',1.5,'Color',pcmap(cond1P,:)})
        hold on
        shadedErrorBar(cond2Time,cond2Traj(:,i),cond2CI(:,i),'lineprops',{'LineWidth',1.5,'Color',pcmap(cond2P,:)})
    end
    
    %Get axis ranges
    VAF = round(trajStruct(1).avgPCA.VAF);
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
       ylabel(['PC ',num2str(plotInd),' (',num2str(VAF(plotInd)),'%)'])
       switch task
           case 'BCI'
               xlabel('time from go cue (ms)')
           case 'planning'
               xlabel('time from target onset (ms)')
           case 'iso'
               xlabel('time from movement onset (ms)')
       end
    end
    
    
    if saveFig
       saveas(gcf,fullfile(saveDir,'UnshiftedPCATimecourses.jpg')); 
    end
    
    %Shifted
    f = figure; f.Position = [100 100 1500 500];
    for i = 1:10
        subplot(2,5,i)
        shadedErrorBar(cond1Time,cond1Traj(:,i),cond1CI(:,i),'lineprops',{'LineWidth',1.5,'Color',pcmap(cond1P,:)})
        hold on
        shadedErrorBar(cond2Time,cond2ShiftTraj(:,i),cond2CI(:,i),'lineprops',{'LineWidth',1.5,'Color',pcmap(cond2P,:)})
    end


    
    %Get axis ranges
    VAF = round(trajStruct(1).avgPCA.VAF);
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
       ylabel(['PC ',num2str(plotInd),' (',num2str(VAF(plotInd)),'%)'])
       switch task
           case 'BCI'
               xlabel('time from go cue (ms)')
           case 'planning'
               xlabel('time from target onset (ms)')
           case 'iso'
               xlabel('time from movement onset (ms)')
       end
    end
        
    if saveFig
       saveas(gcf,fullfile(saveDir,'ShiftedPCATimecourses.jpg')); 
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
     
%% Get optimal alpha
    function alpha = getOptimalAlpha(traj1,traj2,numPts)
        numDims = size(traj1,2);
        alpha = NaN(1,numDims);
        for dim = 1:numDims
            alpha(dim) = (1/numPts)*sum(traj1(1:numPts,dim)-traj2(1:numPts,dim));
        end
    end