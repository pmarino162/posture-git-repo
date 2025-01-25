clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    %saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220906 - compositionality analysis';
    %saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220906 - compositionality analysis\compositionality';
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Figure 8';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%% Construct resultStruct for multiple sessions
    %Earl Reaching
   inputStruct(1).dataset = 'E20210706'; inputStruct(1).task = 'Reach'; inputStruct(1).epoch = 'early';
   inputStruct(2).dataset = 'E20210707'; inputStruct(2).task = 'Reach'; inputStruct(2).epoch = 'early';
   inputStruct(3).dataset = 'E20210708'; inputStruct(3).task = 'Reach'; inputStruct(3).epoch = 'early';
   
   %Nigel Reaching 
    inputStruct(4).dataset = 'N20190222'; inputStruct(4).task = 'Reach'; inputStruct(4).epoch = 'early';
    inputStruct(5).dataset = 'N20190226'; inputStruct(5).task = 'Reach'; inputStruct(5).epoch = 'early';
    inputStruct(6).dataset = 'N20190227'; inputStruct(6).task = 'Reach'; inputStruct(6).epoch = 'early';
    inputStruct(7).dataset = 'N20190228'; inputStruct(7).task = 'Reach'; inputStruct(7).epoch = 'early';

    %Rocky Reaching
    inputStruct(8).dataset = 'R20200221'; inputStruct(8).task = 'Reach'; inputStruct(8).epoch = 'early';
    inputStruct(9).dataset = 'R20200222'; inputStruct(9).task = 'Reach'; inputStruct(9).epoch = 'early';
    
    resultStruct = [];
    for i = 1:numel(inputStruct)
        dataset = inputStruct(i).dataset;
        task = inputStruct(i).task;
        epoch = inputStruct(i).epoch;
        [tempResultStruct] = leaveOneConditionOutCompAnalysis(dataset,task,epoch);
        resultStruct = vertcat(resultStruct,tempResultStruct);
    end
    

%% Create R2 Scatter
    fs = 20;
    for i = 1:numel(resultStruct)
        monkeyList{i} = resultStruct(i).dataset(1);
    end
    figure; hold on
    monkeyInd = 1;
    for monkey = {'E','N','R'}
       tempResultStruct = resultStruct(strcmp(monkeyList,monkey{1,1}));
       R2 = [tempResultStruct.R2];
       meanR2 = mean(R2);
       stdR2 = std(R2);
       numR2 = numel(R2);
       %stdErr = stdR2/sqrt(numR2);
       %errorbar(monkeyInd,meanR2,stdErr);
       plot(monkeyInd,R2,'.','MarkerSize',30,'Color','b');
       monkeyInd = monkeyInd + 1;
    end
    xlim([0 monkeyInd])
    ylim([-.25 1]);
    yticks([0:.25:1]);
    ylabel('R^2')
    xticks([1:3])
    xticklabels({'Monkey E','Monkey N','Monkey R'});
    xtickangle(45)
    ax = gca; 
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
    if saveFig
        saveas(gcf,fullfile(saveDir,'Reach_R2.svg'));
    end
    
    
    %% Visualize results for a dataset
    %Get actual and predicted traj
    tempResultStruct = resultStruct(8);
    actualTrajStruct = tempResultStruct.actualTrajStruct;
    predTrajStruct = tempResultStruct.predTrajStruct;
    
    % Get minimum number of timestamps in condition averages
    minNumTimestamps = size(actualTrajStruct(1).traj,1);  
    
    %Get posture and target lists
    postureList = unique([actualTrajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([actualTrajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(actualTrajStruct(1).traj,2);
    numConditions = size(actualTrajStruct,2);
    
    %Collect both actual and predicted data into X for PCA and marginalization 
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    X = NaN(minNumTimestamps,2*numTargets,numPostures,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).traj;
            predTraj = predTrajStruct([predTrajStruct.posture]==posture & [predTrajStruct.target]==target).traj;
            X(:,targetInd,postureInd,:) = actualTraj(1:minNumTimestamps,:); 
            X(:,targetInd+1,postureInd,:) = predTraj(1:minNumTimestamps,:); 
            targetInd = targetInd + 2;
        end
        postureInd = postureInd+1;
    end
    
    %Get PC space
    [allPCA,score,latent,tsquared,explained,allMu] = pca(reshape(X,[2*numTargets*numPostures*minNumTimestamps,numChannels]));
    
    %Get CI Space
    CIMargTraj = mean(X,[2 3],'omitnan');
    [CIPCA,score,latent,tsquared,explained,mu] = pca(squeeze(CIMargTraj)); 
    
    %Get target space
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj;
    targetMargTraj = squeeze(targetMargTraj);
    targetMargTraj = reshape(targetMargTraj,[2*numTargets*minNumTimestamps,numChannels]);
    [targetPCA,score,latent,tsquared,explained,mu] = pca(targetMargTraj); 
    
    %Get posture space
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1;   
    for posture = postureList
       obsStruct(structInd).label = posture;
       allObs = [];
       tempTrajStruct = actualTrajStruct([actualTrajStruct.posture]==posture);
       for i = 1:size(tempTrajStruct,2)
           traj = tempTrajStruct(i).traj;
           allObs = vertcat(allObs,traj);
%            for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
%                timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
%                traj = tempTrajStruct(i).allSmoothFR(j).traj;
%                if size(traj,1) > minNumTimestamps
%                     traj = traj(1:minNumTimestamps,:);
%                end
%                allObs = vertcat(allObs,traj);
%            end
       end
       obsStruct(structInd).allObs = allObs;
       obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
       structInd = structInd + 1;
    end
    [postureLDA] = doLDA(obsStruct);
    
    %Orthonormalize combinations of axes
    [CPTOrth,~] = qr([postureLDA(:,1),targetPCA(:,1),CIPCA(:,1)]); CPTOrth = CPTOrth(:,1:3);
    %[CPOrth,~] = qr([postureLDA(:,1:2),CIPCA(:,1)]); CPOrth = CPOrth(:,1:3);
    [PTOrth,~] = qr([postureLDA(:,1),targetPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
    [CPTAllOrth,~] = qr([postureLDA(:,1),targetPCA(:,1:2),CIPCA(:,1)]); CPTAllOrth = CPTAllOrth(:,1:4);
        
    %Add Projections to trajStruct
    %totalVar = trace(cov(allTraj));
    for i = 1:size(actualTrajStruct,2)
        %Add average traces to trajStruct
        actualTrajStruct(i).PCA = (actualTrajStruct(i).traj)*allPCA;
        predTrajStruct(i).PCA = (predTrajStruct(i).traj)*allPCA;
        
        actualTrajStruct(i).PTOrth = (actualTrajStruct(i).traj)*PTOrth;
        predTrajStruct(i).PTOrth = (predTrajStruct(i).traj)*PTOrth;
        
%         actualTrajStruct(i).targetPCA = (actualTrajStruct(i).traj)*targetPCA;
%         predTrajStruct(i).targetPCA = (predTrajStruct(i).traj)*targetPCA;
%         
%         actualTrajStruct(i).targetPCA = (actualTrajStruct(i).traj)*targetPCA;
%         predTrajStruct(i).targetPCA = (predTrajStruct(i).traj)*targetPCA;
%         
%         trajStruct(i).PCA.traj = (trajStruct(i).avgSmoothFR.traj-allMu)*allPCA;
%         trajStruct(i).CIPCA.traj = (trajStruct(i).avgSmoothFR.traj-CIMargOffset')*CIPCA;
%         trajStruct(i).targetPCA.traj = trajStruct(i).avgSmoothFR.traj*targetPCA;
%         trajStruct(i).postureLDA.traj = trajStruct(i).avgSmoothFR.traj*postureLDA;
%         trajStruct(i).CPTOrth.traj = trajStruct(i).avgSmoothFR.traj*CPTOrth;
%         %trajStruct(i).CPOrth.traj = trajStruct(i).avgSmoothFR.traj*CPOrth;
%         trajStruct(i).PTOrth.traj = trajStruct(i).avgSmoothFR.traj*PTOrth;
%         trajStruct(i).CPTAllOrth.traj = trajStruct(i).avgSmoothFR.traj*CPTAllOrth;
%         %Get VAF
%         trajStruct(i).PCA.VAF =  100.*(diag(cov(allTraj*allPCA))')./totalVar;
%         trajStruct(i).PTOrth.VAF =  100.*(diag(cov(allTraj*PTOrth))')./totalVar;
%         trajStruct(i).postureLDA.VAF = 100.*(diag(cov(allTraj*postureLDA))')./totalVar;
    end
        
     
    

    
    
    %3d Posture-Target views of actual and predicted data
    for i = 1:2
       if i == 1
            tempTrajStruct = actualTrajStruct;
       else
            tempTrajStruct = predTrajStruct;
       end
       figure; hold on
       xDim = 2; yDim = 3; zDim = 1;
        for posture = postureList
            for target = targetList
                if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
                    traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth; 
                    if posture <= 5
                        plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    else
                        plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(mod(posture,5),:),'LineWidth',2);
                        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(mod(posture,5),:),'MarkerFaceColor',pcmap(mod(posture,5),:));
                        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(mod(posture,5),:),'MarkerFaceColor',pcmap(mod(posture,5),:));
                    end
                end
            end
       end
       xlabel('Posture LDA')
       ylabel('Target PC1')
       zlabel('Target PC2')
       view([10 20])
       xticks([]); yticks([]); zticks([]);
       grid on
        ax = gca; 
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
       if saveFig
           if i==1
                saveas(gcf,fullfile(saveDir,[dataset,'_ActualPTOrth.svg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_ActualPTOrth.fig']));
           elseif i==2
                saveas(gcf,fullfile(saveDir,[dataset,'_ModelPTOrth.svg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_ModelPTOrth.fig']));
           end
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