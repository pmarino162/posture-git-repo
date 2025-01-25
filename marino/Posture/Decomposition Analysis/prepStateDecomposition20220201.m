clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    task = 'planning';
    switch task
        case 'planning'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\Planning Figures';
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
        case 'planning'
        [Data] = loadEarlData20210706_20211210; 
        [Data] = subselectForMarginalization(Data,task);
         trialInclStates(1).trialName = {'GridReaching'};
        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        trialInclStates(1).inclStates = {{'state','Delay','first',0},{'state','Target Acquire','first',0}};
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
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
      
%% Create Prep State obs struct
    obsStruct = struct('posture',[],'target',[],'numObs',[],'allObs',[],'grandMean',[]);
    structInd = 1;   
    for posture = postureList
       for target = targetList
           obsStruct(structInd).posture = posture;
           obsStruct(structInd).target = target;
           allObs = [];
           tempTrajStruct = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target);
           allSmoothFR = tempTrajStruct(1).allSmoothFR;
           for i = 1:size(allSmoothFR,2)
               timestamps = allSmoothFR(i).timestamps;
               traj = allSmoothFR(i).traj;
               traj = mean(traj(11:end,:));
               allObs = vertcat(allObs,traj);
           end
           obsStruct(structInd).allObs = allObs;
           obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
           obsStruct(structInd).grandMean = mean(allObs);
           structInd = structInd + 1;
       end
    end
    
%% Do PCA on condition averages for visualization; get CIs and Var expl
%     allAvgs = NaN(minNumTimestamps*size(trajStruct,2),size(trajStruct(1).allSmoothFR(1).traj,2));
%     j = 1;
%     for i = 1:size(trajStruct,2)
%         allAvgs(j:j+minNumTimestamps-1,:) = trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps,:);
%         j = j+minNumTimestamps;
%     end
    allAvgs = vertcat(obsStruct.allObs);
    [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);       
    for i = 1:size(obsStruct,2)
        obsStruct(i).allObsPC = (obsStruct(i).allObs-mu)*coeff(:,1:10);
        obsStruct(i).grandMeanPC = mean(obsStruct(i).allObsPC);
    end
    
%% Visualize
    for i = 1:size(obsStruct,2)
       posture = obsStruct(i).posture;
       target = obsStruct(i).target;
       allObs = obsStruct(i).allObs;
       allObs = (allObs-mu)*coeff;
       plot3(allObs(:,1),allObs(:,2),allObs(:,3),'.','MarkerSize',20,'Color',tcmap(target,:));
       hold on;
    end
    
%% Peform LDA by posture on all data 
    postureObsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;
    for posture = postureList
        tempObsStruct = obsStruct([obsStruct.posture]==posture);
        postureObsStruct(structInd).label = posture;
        postureObsStruct(structInd).numObs = sum([tempObsStruct.numObs]);
        postureObsStruct(structInd).allObs = vertcat(tempObsStruct.allObs);
        structInd = structInd + 1;
    end
    [postureLDA] = doLDA(postureObsStruct);

%% Peform LDA by target on all data 
    targetObsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;
    for target = targetList
        tempObsStruct = obsStruct([obsStruct.target]==target);
        targetObsStruct(structInd).label = target;
        targetObsStruct(structInd).numObs = sum([tempObsStruct.numObs]);
        targetObsStruct(structInd).allObs = vertcat(tempObsStruct.allObs);
        structInd = structInd + 1;
    end
    [targetLDA] = doLDA(targetObsStruct);
    
    p2TargetObsStruct = struct('label',[],'numObs',[],'allObs',[]);
    p2ObsStruct = obsStruct([obsStruct.posture]==2 | [obsStruct.posture]==3);
    structInd = 1;
    for target = targetList
        tempObsStruct = p2ObsStruct([p2ObsStruct.target]==target);
        p2TargetObsStruct(structInd).label = target;
        p2TargetObsStruct(structInd).numObs = sum([tempObsStruct.numObs]);
        p2TargetObsStruct(structInd).allObs = vertcat(tempObsStruct.allObs);
        structInd = structInd + 1;
    end
    [p2targetLDA] = doLDA(p2TargetObsStruct);
    
%% Do LDA by target on PC proj
    p2PCTargetObsStruct = struct('label',[],'numObs',[],'allObs',[]);
    p2ObsStruct = obsStruct([obsStruct.posture]==1);
    structInd = 1;
    for target = targetList
        tempObsStruct = p2ObsStruct([p2ObsStruct.target]==target);
        p2PCTargetObsStruct(structInd).label = target;
        p2PCTargetObsStruct(structInd).numObs = sum([tempObsStruct.numObs]);
        p2PCTargetObsStruct(structInd).allObs = vertcat(tempObsStruct.allObsPC);
        structInd = structInd + 1;
    end
    [p2PCtargetLDA] = doLDA(p2PCTargetObsStruct);

%% Do rings plot a la adam
    [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:3);
    %Grand means in p2 target/posture space
    figure
    for posture = postureList
       tempObsStruct = obsStruct([obsStruct.posture]==posture);
       allGrandMeans = vertcat(tempObsStruct.grandMean);
       allGrandMeans = vertcat(allGrandMeans,allGrandMeans(1,:));
       for i = 1:size(tempObsStruct,2)
           target = tempObsStruct(i).target;
           grandMean = tempObsStruct(i).grandMean;
           LDAproj = grandMean*postTargOrth;
           plot3(LDAproj(:,1),LDAproj(:,2),LDAproj(:,3),'o','MarkerSize',10,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',pcmap(posture,:));
           hold on;
       end
       LDAproj = allGrandMeans*postTargOrth;
       plot3(LDAproj(:,1),LDAproj(:,2),LDAproj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
       xlabel('1');ylabel('2');zlabel('3')
    end

%% Check to see if things hold up using only one posture (PCs)
    %Try w PCS
    [targOrth,~] = qr(p2PCtargetLDA(:,1:2));
    targOrth = targOrth(:,1:2);
    %Grand means in p2 target space
    figure
    for posture = postureList
       tempObsStruct = obsStruct([obsStruct.posture]==posture);
       allGrandMeans = vertcat(tempObsStruct.grandMeanPC);
       allGrandMeans = vertcat(allGrandMeans,allGrandMeans(1,:));
       for i = 1:size(tempObsStruct,2)
           target = tempObsStruct(i).target;
           grandMean = tempObsStruct(i).grandMeanPC;
           LDAproj = grandMean*targOrth;
           plot(LDAproj(:,1),LDAproj(:,2),'o','MarkerSize',10,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',pcmap(posture,:));
           hold on;
       end
       LDAproj = allGrandMeans*targOrth;
       plot(LDAproj(:,1),LDAproj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
       xlabel('1');ylabel('2');zlabel('3')
    end

%% Plot top PCs and Posture 1
    figure
    for i = 1:size(obsStruct,2)
       posture = obsStruct(i).posture;
       target = obsStruct(i).target;
       allObs = obsStruct(i).allObs;
       PCs = (allObs-mu)*coeff;
       pLDA = allObs*postureLDA;
       plot3(PCs(:,1),PCs(:,2),pLDA(:,1),'o','MarkerSize',10,'MarkerFaceColor',pcmap(posture,:),'MarkerEdgeColor',tcmap(target,:));
       hold on;
       xlabel('1');ylabel('2');zlabel('3')
    end
    
    figure
    for i = 1:size(obsStruct,2)
       posture = obsStruct(i).posture;
       target = obsStruct(i).target;
       allObs = obsStruct(i).allObs;
       PCs = (allObs-mu)*coeff;
       pLDA = allObs*postureLDA;
       plot3(PCs(:,1),PCs(:,2),pLDA(:,1),'o','MarkerSize',10,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',pcmap(posture,:));
       hold on;
       xlabel('1');ylabel('2');zlabel('3')
    end
    
%% Posture and target
    figure
    for i = 1:size(obsStruct,2)
       posture = obsStruct(i).posture;
       target = obsStruct(i).target;
       allObs = obsStruct(i).allObs;
       tLDA = allObs*targetLDA;
       pLDA = allObs*postureLDA;
       plot3(tLDA(:,1),tLDA(:,2),pLDA(:,1),'o','MarkerSize',10,'MarkerFaceColor',pcmap(posture,:),'MarkerEdgeColor',tcmap(target,:));
       hold on;
       xlabel('1');ylabel('2');zlabel('3')
    end
    
    figure
    for i = 1:size(obsStruct,2)
       posture = obsStruct(i).posture;
       target = obsStruct(i).target;
       allObs = obsStruct(i).allObs;
       tLDA = allObs*targetLDA
       pLDA = allObs*postureLDA;
       plot3(tLDA(:,1),tLDA(:,2),pLDA(:,1),'o','MarkerSize',10,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',pcmap(posture,:));
       hold on;
       xlabel('1');ylabel('2');zlabel('3')
    end

%% Look at all data, grand means, in p2 and full LDA space
%     [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
%     postTargOrth = postTargOrth(:,1:3);
% 
%     %Grand means in posture/ target space
%     figure
%     for i = 1:size(obsStruct,2)
%        posture = obsStruct(i).posture;
%        target = obsStruct(i).target;
%        grandMean = obsStruct(i).grandMean;
%        LDAproj = grandMean*postTargOrth;
%        plot3(LDAproj(:,1),LDAproj(:,2),LDAproj(:,3),'o','MarkerSize',10,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',pcmap(posture,:));
%        hold on;
%        xlabel('1');ylabel('2');zlabel('3')
%     end

    [postTargOrth,~] = qr([postureLDA(:,1),p2targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:3);
    %Grand means in p2 target/posture space
    figure
    for posture = postureList
       tempObsStruct = obsStruct([obsStruct.posture]==posture);
       allGrandMeans = vertcat(tempObsStruct.grandMean);
       allGrandMeans = vertcat(allGrandMeans,allGrandMeans(1,:));
       for i = 1:size(tempObsStruct,2)
           target = tempObsStruct(i).target;
           grandMean = tempObsStruct(i).grandMean;
           LDAproj = grandMean*postTargOrth;
           plot3(LDAproj(:,1),LDAproj(:,2),LDAproj(:,3),'o','MarkerSize',10,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',pcmap(posture,:));
           hold on;
       end
       LDAproj = allGrandMeans*postTargOrth;
       plot3(LDAproj(:,1),LDAproj(:,2),LDAproj(:,3),'Color',pcmap(posture,:),'LineWidth',2);
       xlabel('1');ylabel('2');zlabel('3')
    end
    
    [targOrth,~] = qr(p2targetLDA(:,1:2));
    targOrth = targOrth(:,1:3);
    %Grand means in p2 target space
    figure
    for posture = postureList
       tempObsStruct = obsStruct([obsStruct.posture]==posture);
       allGrandMeans = vertcat(tempObsStruct.grandMean);
       allGrandMeans = vertcat(allGrandMeans,allGrandMeans(1,:));
       for i = 1:size(tempObsStruct,2)
           target = tempObsStruct(i).target;
           grandMean = tempObsStruct(i).grandMean;
           LDAproj = grandMean*targOrth;
           plot(LDAproj(:,1),LDAproj(:,2),'o','MarkerSize',10,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',pcmap(posture,:));
           hold on;
       end
       LDAproj = allGrandMeans*targOrth;
       plot(LDAproj(:,1),LDAproj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
       xlabel('1');ylabel('2');zlabel('3')
    end
    
    %Try w PCS
    [targOrth,~] = qr(p2PCtargetLDA(:,1:2));
    targOrth = targOrth(:,1:2);
    %Grand means in p2 target space
    figure
    for posture = postureList
       tempObsStruct = obsStruct([obsStruct.posture]==posture);
       allGrandMeans = vertcat(tempObsStruct.grandMeanPC);
       allGrandMeans = vertcat(allGrandMeans,allGrandMeans(1,:));
       for i = 1:size(tempObsStruct,2)
           target = tempObsStruct(i).target;
           grandMean = tempObsStruct(i).grandMeanPC;
           LDAproj = grandMean*targOrth;
           plot(LDAproj(:,1),LDAproj(:,2),'o','MarkerSize',10,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',pcmap(posture,:));
           hold on;
       end
       LDAproj = allGrandMeans*targOrth;
       plot(LDAproj(:,1),LDAproj(:,2),'Color',pcmap(posture,:),'LineWidth',2);
       xlabel('1');ylabel('2');zlabel('3')
    end
    
    
%     %allData in p2 target space
%     figure
%     for i = 1:size(obsStruct,2)
%        posture = obsStruct(i).posture;
%        target = obsStruct(i).target;
%        allObs = obsStruct(i).allObs;
%        LDAproj = allObs*postTargOrth;
%        plot3(LDAproj(:,1),LDAproj(:,2),LDAproj(:,3),'o','MarkerSize',10,'MarkerFaceColor',tcmap(target,:),'MarkerEdgeColor',pcmap(posture,:));
%        hold on;
%        xlabel('1');ylabel('2');zlabel('3')
%     end
    
%% Peform LDA by target on all data - around 250ms after cue onset
    %Get number of points for each condition
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;
    for target = targetList
       obsStruct(structInd).label = target;
       allObs = [];
       switch task
           case 'BCI'
                tempTrajStruct = trajStruct([trajStruct.target]==target);
                timeToUse = 125;
                numOnEitherSide = 1;
           case 'iso'
                tempTrajStruct = trajStruct([trajStruct.target]==target);
                timeToUse = -150;
                numOnEitherSide = 0;
           case 'reaching'
                tempTrajStruct = trajStruct([trajStruct.target]==target);
                timeToUse = 250;
                numOnEitherSide = 0;
           case 'planning'
                tempTrajStruct = trajStruct([trajStruct.target]==target);
                timeToUse = 150;
                numOnEitherSide = 0;
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
     