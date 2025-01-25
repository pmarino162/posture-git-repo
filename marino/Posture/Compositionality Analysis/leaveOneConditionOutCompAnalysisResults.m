clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = true;
    %saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220906 - compositionality analysis';
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220906 - compositionality analysis\compositionality';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%% Construct resultStruct for multiple sessions
%     %Earl BCI
%     inputStruct(1).dataset = 'E20200316'; inputStruct(1).task = 'BCI'; inputStruct(1).epoch = 'all';
%     inputStruct(2).dataset = 'E20200317'; inputStruct(2).task = 'BCI'; inputStruct(2).epoch = 'all';
%     inputStruct(3).dataset = 'E20200318'; inputStruct(3).task = 'BCI'; inputStruct(3).epoch = 'all';

%      %Earl Reaching
%    inputStruct(1).dataset = 'E20210706'; inputStruct(1).task = 'Reach'; inputStruct(1).epoch = 'all';
%    inputStruct(2).dataset = 'E20210707'; inputStruct(2).task = 'Reach'; inputStruct(2).epoch = 'all';
%     inputStruct(3).dataset = 'E20210708'; inputStruct(3).task = 'Reach'; inputStruct(3).epoch = 'all';
%     %Earl Multiple Tasks
%     inputStruct(1).dataset = 'E20200314'; inputStruct(1).task = 'BCI'; inputStruct(1).epoch = 'all';
%     inputStruct(2).dataset = 'E20200314'; inputStruct(2).task = 'Reach'; inputStruct(2).epoch = 'all';
%     inputStruct(3).dataset = 'E20200314'; inputStruct(3).task = 'Iso'; inputStruct(3).epoch = 'all';


%      %Nigel BCI
%      inputStruct(1).dataset = 'N20171215'; inputStruct(1).task = 'BCI'; inputStruct(1).epoch = 'all';
%      inputStruct(2).dataset = 'N20180221'; inputStruct(2).task = 'BCI'; inputStruct(2).epoch = 'all';

%     %Nigel Reaching 
%     inputStruct(1).dataset = 'N20190222'; inputStruct(1).task = 'Reach'; inputStruct(1).epoch = 'all';
%     inputStruct(2).dataset = 'N20190226'; inputStruct(2).task = 'Reach'; inputStruct(2).epoch = 'all';
%     inputStruct(3).dataset = 'N20190227'; inputStruct(3).task = 'Reach'; inputStruct(3).epoch = 'all';
%     inputStruct(4).dataset = 'N20190228'; inputStruct(4).task = 'Reach'; inputStruct(4).epoch = 'all';

%     %Rocky Reaching
%     inputStruct(1).dataset = 'R20200221'; inputStruct(1).task = 'Reach'; inputStruct(1).epoch = 'all';
%     inputStruct(2).dataset = 'R20200222'; inputStruct(2).task = 'Reach'; inputStruct(2).epoch = 'all';


    resultStruct = [];
    for i = 1:numel(inputStruct)
        dataset = inputStruct(i).dataset;
        task = inputStruct(i).task;
        epoch = inputStruct(i).epoch;
        [tempResultStruct] = leaveOneConditionOutCompAnalysis(dataset,task,epoch);
        resultStruct = vertcat(resultStruct,tempResultStruct);
    end
    
%% Get monkey and task info for plotting
    monkey = resultStruct(1).dataset(1);
    task = resultStruct(1).task;
    
%% Create R2 Histogram 
figure ; hold on
    histogram([resultStruct.CIR2],'BinWidth',.1);
    histogram([resultStruct.CIPR2],'BinWidth',.1);
    histogram([resultStruct.R2],'BinWidth',.1);
    
%     histogram([resultStruct.CIFullR2],'BinWidth',.1);
%     histogram([resultStruct.CIPFullR2],'BinWidth',.1);
%     histogram([resultStruct.fullR2],'BinWidth',.1);
%     
    xlim([0 1])
    xlabel('R^2')
    ylabel('Num Sessions')
    legend({'CI','CI+P','CI+P+T'})
    if saveFig
        saveas(gcf,fullfile(saveDir,[monkey,'_',task,'_R2Hist.svg']));
    end
    
%% Create mean dist histogram 
figure ; hold on
%     histogram([resultStruct.noModelDist],'BinWidth',.1);
%     histogram([resultStruct.CIFullDist],'BinWidth',.1);
%     histogram([resultStruct.CIPFullDist],'BinWidth',.1);
%     histogram([resultStruct.fullFullDist],'BinWidth',.1);
    
    histogram([resultStruct.noModelDist],'BinWidth',.1);
    histogram([resultStruct.CIDist],'BinWidth',.1);
    histogram([resultStruct.CIPDist],'BinWidth',.1);
    histogram([resultStruct.fullDist],'BinWidth',.1);
%     
    %xlim([0 5])
    ax = gca;
    xlims = ax.XLim;
    xlim([0,xlims(2)])
    xlabel('Mean Dist (std)')
    ylabel('Num Sessions')
    legend({'No Model','CI','CI+P','CI+P+T'})
    
    if saveFig
        saveas(gcf,fullfile(saveDir,[monkey,'_',task,'_meanDistHist.svg']));
    end
%% Willett Pie 
    tempResultStruct = resultStruct(1);
    figure
    CIVarPct = tempResultStruct.CIVarPct;
    postureVarPct = tempResultStruct.postureVarPct;
    targetVarPct = tempResultStruct.targetVarPct;
    intVarPct = tempResultStruct.intVarPct;
    labels = {['Condition Invariant ',num2str(round(CIVarPct)),'%'],['Target ',num2str(round(targetVarPct)),'%'],['Posture ',num2str(round(postureVarPct)),'%'],['PxT Interaction ',num2str(round(intVarPct)),'%']};
    pie([CIVarPct,targetVarPct,postureVarPct,intVarPct],[1 1 1 1],labels)
    if saveFig
        saveas(gcf,fullfile(saveDir,[monkey,'_',task,'_WillettPie.svg']));
    end
    
    %% Visualize results for a dataset
    %Get actual and predicted traj
    tempResultStruct = resultStruct(1);
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
        
     
    %3d PCA views of actual and predicted data
    for i = 1:2
       if i == 1
            tempTrajStruct = actualTrajStruct;
       else
            tempTrajStruct = predTrajStruct;
       end
       figure; hold on
       xDim = 1; yDim = 2; zDim = 3;
        for posture = [postureList(1),postureList(end)]
            for target = targetList
                if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
                    traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PCA; 
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
       xlabel('PC 1')
       ylabel('PC 2')
       zlabel('PC 3')
       if i == 1
           title('Actual')
       else
           title('Model')
       end
       view([10 20])
       
       if saveFig
           if i==1
                saveas(gcf,fullfile(saveDir,[dataset,'_ActualPC.svg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_ActualPC.fig']));
           elseif i==2
                saveas(gcf,fullfile(saveDir,[dataset,'_ModelPC.svg']));
                saveas(gcf,fullfile(saveDir,[dataset,'_ModelPC.fig']));
           end
       end
       
    end
    

    
    
    %3d Posture-Target views of actual and predicted data
    for i = 1:2
       if i == 1
            tempTrajStruct = actualTrajStruct;
       else
            tempTrajStruct = predTrajStruct;
       end
       figure; hold on
       xDim = 1; yDim = 2; zDim = 3;
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
       if i == 1
           title('Actual')
       else
           title('Model')
       end
       
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
    
%     %
%     for i = 1:4
%         if i == 1 || i == 3
%             tempTrajStruct = actualTrajStruct;
%         else
%             tempTrajStruct = predTrajStruct;
%         end
%         figure; hold on
%         xDim = 1; yDim = 2; zDim = 3;
%         for posture = postureList
%             for target = targetList
%                 if any([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)
%                     if i == 1 || i == 2
%                         traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PCA; 
%                     else
%                         traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).PTOrth; 
%                     end
%                     plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
%                     plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
%                     plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
%                 end
%             end
%         end
%     end
    
    %3d Posture-Target view of predicted data
    
    
    
%     
%     %Plot in top PCs
%     posture = 1;
%     figure
%     for target = 1:8
%         actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).traj;
%         actualTraj = actualTraj*coeff(:,1:10);
%         for PC = 1:10
%             subplot(2,5,PC)
%             plot(actualTraj(:,PC),'Color',tcmap(target,:));
%             hold on
%         end
%     end
%     
%     figure
%     posture = 1;
%     for target = 1:8
%         predTraj = predTrajStruct([predTrajStruct.posture]==posture & [predTrajStruct.target]==target).traj;
%         predTraj = predTraj*coeff(:,1:10);
%         for PC = 1:10
%             subplot(2,5,PC)
%             plot(predTraj(:,PC),'Color',tcmap(target,:));
%             hold on
%         end
%     end
%     
%     
%     
%     
%     figure
%     target = 5;
%     for posture = 1:2
%         actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).traj;
%         actualTraj = actualTraj*coeff(:,1:10);
%         for PC = 1:10
%             subplot(2,5,PC)
%             plot(actualTraj(:,PC),'Color',pcmap(posture,:));
%             hold on
%         end
%     end
%     
%     
%     figure
%     target = 1;
%     for posture = 1:2
%         predTraj = predTrajStruct([predTrajStruct.posture]==posture & [predTrajStruct.target]==target).traj;
%         predTraj = predTraj*coeff(:,1:10);
%         for PC = 1:10
%             subplot(2,5,PC)
%             plot(predTraj(:,PC),'Color',pcmap(posture,:));
%             hold on
%         end
%     end
%     
%     figure;
%     target = 2; posture = 2;
%     actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).traj;
%     actualTraj = actualTraj*coeff(:,1:10);
%     predTraj = predTrajStruct([predTrajStruct.posture]==posture & [predTrajStruct.target]==target).traj;
%     predTraj = predTraj*coeff(:,1:10);
%     for PC = 1:10
%         subplot(2,5,PC); hold on
%         plot(actualTraj(:,PC),'Color','b');
%         plot(predTraj(:,PC),'Color','r');
%     end
%     
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