clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20221123\Figures\Figure 2\Alternate Upper Panel Projections 20221219';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
        
%% Run loop for each dataset    
    for datasetList = {'E20210706'}%{'R20201020'}%{'E20200116' 'E20210706','N20180221'} 
        %datasetList = {'E20200317','N20180221'} 
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            case 'E20200317'
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            case 'E20200116'
                trialInclStates(1).trialName = {'IsometricForce_1D'};   
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',250}};
            case 'E20210706'
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
        end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        %Keep postures 1-4 only for earl
        switch dataset
            case {'E20210706'}
                trajStruct = trajStruct(ismember([trajStruct.posture],[1,2,3,4]));
        end
        
        % For each timecourse in trajstruct, get single point
        for i = 1:size(trajStruct,2)
            %Individual trials
            for j = 1:size(trajStruct(i).allSmoothFR,2)
               trajStruct(i).allSmoothFR(j).trialAvg = mean(trajStruct(i).allSmoothFR(j).traj); 
            end
        end
    
        % Get minimum number of timestamps in condition averages
        numTimestamps = [];
        for i = 1:size(trajStruct,2)
            numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
        end
        [minNumTimestamps,i] = min(numTimestamps);
        minNumTimestamps = 10;
        
        
        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);
        
        %Get obsStruct
        obsStruct = struct('label',[],'numObs',[],'allObs',[]);
        label = 1;
        for posture = postureList
            allObs = [];
            for target = targetList
                allSmoothFR = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allSmoothFR;
                for i = 1:size(allSmoothFR,2)
                   obs = allSmoothFR(i).traj;
                   allObs = vertcat(obs,allObs);
                end
            end
            obsStruct(label).label = label;
            obsStruct(label).allObs = allObs;
            obsStruct(label).numObs = size(allObs,1);
            label = label + 1;
        end

        %Do PCA on observations
        allObs = vertcat(obsStruct.allObs);
        [coeff,score,latent,tsquared,explained,mu] = pca(allObs);
    
        %Replace obsStruct observations w PC projections
        numPCs = 20;
        for i = 1:numel(obsStruct)
           obsStruct(i).allObs = obsStruct(i).allObs*coeff(:,1:numPCs);       
        end

        %Do LDA on obsStruct
        [postureLDA] = doLDA(obsStruct);

        %Combine PCA and LDA to get posture Axis 
        postureLDA = coeff(:,1:numPCs)*postureLDA;
        
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
        postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
    
        %Do PCA on Marginalizations    
        CIMargTraj = squeeze(CIMargTraj);
            CIMargOffset = squeeze(CIMargOffset);
            [CIPCA,score,latent,tsquared,explained,mu] = pca(CIMargTraj); 
        targetMargTraj = squeeze(targetMargTraj);
            targetMargTraj = reshape(targetMargTraj,[numTargets*minNumTimestamps,numChannels]);
            [targetPCA,score,latent,tsquared,explained,mu] = pca(targetMargTraj); 
        %postureMargOffset = squeeze(postureMargOffset);
        postureMargTraj = squeeze(postureMargTraj);
            postureMargTraj = reshape(postureMargTraj,[numPostures*minNumTimestamps,numChannels]);
            [posturePCA,score,latent,tsquared,explained,mu] = pca(postureMargTraj); 

        %Define orthonormalized spaces
        
        postureLDA(:,1) = -1.*postureLDA(:,1);
        posturePCA(:,1) = -1.*posturePCA(:,1);
        [PTOrthPCA,~] = qr([posturePCA(:,1),targetPCA(:,1:2)]); PTOrthPCA = PTOrthPCA(:,1:3);
        [PTOrthLDAPCA,~] = qr([postureLDA(:,1),targetPCA(:,1:2)]); PTOrthLDAPCA = PTOrthLDAPCA(:,1:3);
            
        %Add Projections to trajStruct
        totalVar = trace(cov(allTraj));
        for i = 1:size(trajStruct,2)
            %Add average traces to trajStruct
            trajStruct(i).PCA.traj = (trajStruct(i).avgSmoothFR.traj-allMu)*allPCA;
            trajStruct(i).PTOrthLDAPCA.traj = trajStruct(i).avgSmoothFR.traj*PTOrthLDAPCA;
            trajStruct(i).PTOrthPCA.traj = trajStruct(i).avgSmoothFR.traj*PTOrthPCA;
            %Get VAF
            trajStruct(i).PCA.VAF =  100.*(diag(cov(allTraj*allPCA))')./totalVar;
            trajStruct(i).PTOrthLDAPCA.VAF =  100.*(diag(cov(allTraj*PTOrthLDAPCA))')./totalVar;
            trajStruct(i).PTOrthPCA.VAF =  100.*(diag(cov(allTraj*PTOrthPCA))')./totalVar;
        end
        
        %Plot Parameters
        fs = 14;
        
        % PT-Orth; scaled & new axis labels
        ICs = zeros(1,3); ICind = 1;
        figure; hold on
        xDim = 2; yDim = 3; zDim = 1;
        for posture = postureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(1:minNumTimestamps,:);
                    traj = traj*PTOrthLDAPCA;
                    ICs(ICind,:) = [traj(1,xDim),traj(1,yDim),traj(1,zDim)];
                    plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    ICind = ICind + 1;
                end
            end
        end
        ax = gca; xticklabels({}); yticklabels({}); zticklabels({}); set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        VAF = round(trajStruct(1).PTOrthLDAPCA.VAF);
        xlabel(['Dim 1 (',num2str(VAF(xDim)),'%)']); ylabel(['Dim 2 (',num2str(VAF(yDim)),'%)']); zlabel(['Dim 3 (',num2str(VAF(zDim)),'%)']);
        
        % PT-Orth; unscaled & new axis labels
        axis equal
        grid on; 
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'PTOrthLDAPCA_unscaled.svg']));
        end
        
        % Add z quiver
        zLen = ax.ZLim(2)-ax.ZLim(1);
        meanIC = mean(ICs,1);
        offset = zLen*.5;
        q = quiver3(meanIC(1),meanIC(2),meanIC(3)+offset,0,0,-1,5,'^','MarkerFaceColor','k');
        q.Color = [0 0 0]; q.LineWidth = 3;
        grid on; 
        %Add x and y quiver
        xLen = ax.XLim(2)-ax.XLim(1);
        offset = xLen*.5;
        q = quiver3(meanIC(1)+offset,meanIC(2),meanIC(3),-1,0,0,5,'>','MarkerFaceColor','k');
        q.Color = [0 0 0]; q.LineWidth = 3;
        
        yLen = ax.YLim(2)-ax.YLim(1);
        offset = yLen*.5;
        q = quiver3(meanIC(1),meanIC(2)-offset,meanIC(3),0,1,0,5,'>','MarkerFaceColor','k');
        q.Color = [0 0 0]; q.LineWidth = 3;
        grid on;        
        
        
        % PT-Orth; scaled & new axis labels
        ICs = zeros(1,3); ICind = 1;
        figure; hold on
        xDim = 2; yDim = 3; zDim = 1;
        for posture = postureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(1:minNumTimestamps,:);
                    traj = traj-CIMargTraj(1:minNumTimestamps,:);
                    traj = traj*PTOrthPCA;
                    ICs(ICind,:) = [traj(1,xDim),traj(1,yDim),traj(1,zDim)];
                    plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    ICind = ICind + 1;
                end
            end
        end
        grid on; ax = gca; xticklabels({}); yticklabels({}); zticklabels({}); set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        VAF = round(trajStruct(1).PTOrthPCA.VAF);
        xlabel(['Dim 1 (',num2str(VAF(xDim)),'%)']); ylabel(['Dim 2 (',num2str(VAF(yDim)),'%)']); zlabel(['Dim 3 (',num2str(VAF(zDim)),'%)']);
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'isoForceTraj.svg']));
        end
        
        
        % PT-Orth; unscaled & new axis labels
        axis equal
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'PTOrthLDAPCA_unscaled.svg']));
        end
        % Add z quiver
        zLen = ax.ZLim(2)-ax.ZLim(1);
        meanIC = mean(ICs,1);
        offset = zLen*.5;
        q = quiver3(meanIC(1),meanIC(2),meanIC(3)+offset,0,0,-1,5,'^','MarkerFaceColor','k');
        q.Color = [0 0 0]; q.LineWidth = 3;
        
        %Add x and y quiver
        xLen = ax.XLim(2)-ax.XLim(1);
        offset = xLen*.5;
        q = quiver3(meanIC(1)+offset,meanIC(2),meanIC(3),-1,0,0,5,'>','MarkerFaceColor','k');
        q.Color = [0 0 0]; q.LineWidth = 3;
        
        yLen = ax.YLim(2)-ax.YLim(1);
        offset = yLen*.5;
        q = quiver3(meanIC(1),meanIC(2)-offset,meanIC(3),0,1,0,5,'>','MarkerFaceColor','k');
        q.Color = [0 0 0]; q.LineWidth = 3;
        
        % Top 3 PCs orth
        switch dataset
            case {'E20200317','E20200116'}
                plotTargetList = [1,5];
                plotPostureList = [1,5];
            case 'E20210706'
                plotPostureList = [1,4];
        end
        figure; hold on
        xDim = 1; yDim = 2; zDim = 3;
        for posture = plotPostureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.traj(1:minNumTimestamps,:); 
                    plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                end
            end
        end           
        grid on; xticklabels({}); yticklabels({}); zticklabels({}); 
        VAF = round(trajStruct(1).PCA.VAF);
        xlabel(['PC 1 (',num2str(VAF(xDim)),'%)']); ylabel(['PC 2 (',num2str(VAF(yDim)),'%)']); zlabel(['PC 3 (',num2str(VAF(zDim)),'%)'])
        view([20 10])
        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PCA_3d.svg']));
        end
        
        % Top PC time courses
        switch dataset
            case {'E20200317','E20210706'}
                plotTargetList = [1,5];
                plotPostureList = [1,4];
            case 'E20200116'
                plotTargetList = [3,7];
                
        end
        f = figure; f.Position = [200 200 300 500];
        for posture = plotPostureList
            for target = plotTargetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PCA.traj(1:minNumTimestamps,:); 
                    time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps(1,1:minNumTimestamps); 
                    for dim = 1:10
                       subplot(5,2,dim) 
                       plot(time,traj(:,dim),'Color',pcmap(posture,:),'LineWidth',2);
                       hold on;
                       xticklabels({}); yticklabels({});
                    end
                end
            end
        end    
        
       linkaxes
        
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PCA_timecourses.svg']));
        end
        
    end
    
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
         