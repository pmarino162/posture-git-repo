clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20230203\orthographic dpca pop traj';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Get trajStruct
    %Load data
    dataset = 'E20200316';
    [Data,zScoreParams] = loadData(dataset);
    
    %Get trajStruct
    binWidth = 25;
    kernelStdDev = 25;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
    
    %Execution
    trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',200}};
    exTrajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    
    %Touch bar hold
    trialInclStates(1).inclStates = {{'state','Center Target','first',-150},{'state','Center Target','first',50}};
    tbHoldTrajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    
%% do dPCA
    %Manually enter number of points to consider (execution)
   	numPts = 8;
    
    % Get minimum number of timestamps in condition averages (tb Hold)
    numHoldPts = [];
    for i = 1:size(tbHoldTrajStruct,2)
        numHoldPts(i) = length(tbHoldTrajStruct(i).avgSmoothFR.timestamps);
    end
    [numHoldPts,~] = min(numHoldPts);
        
    %Get posture and target lists
    postureList = unique([exTrajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([exTrajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(exTrajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(exTrajStruct,2);

    %Form X, containing trial-averaged data for each condition
    X = NaN(numPts,numTargets,numPostures,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.traj;
            X(:,targetInd,postureInd,:) = traj(1:numPts,:); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
        
    %Do PCA on all condition averages
    allTraj = reshape(X,[numTargets*numPostures*numPts,numChannels]);
    mu = mean(allTraj,1);
    %[allPCA,score,latent,tsquared,explained,allMu] = pca(allTraj);
        
    Xdpca = permute(X,[4,3,2,1]);

    combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
    margNames = {'P', 'T', 'CI', 'PTI'};
    margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

    N = numChannels;
    P = numPostures;
    T = numTargets;
    timePts = numPts;
    time = exTrajStruct(1).avgSmoothFR.timestamps(1:numPts);
    timeEvents = [];

    [W,V,whichMarg] = dpca(Xdpca, 10, ...
        'combinedParams', combinedParams);

    explVar = dpca_explainedVariance(Xdpca, W, V, ...
        'combinedParams', combinedParams);

    % Add Projections to trajStruct; get CI
    totalVar = trace(cov(allTraj));
    for i = 1:size(exTrajStruct,2)
        %Add average traces to trajStruct
        exTrajStruct(i).avgdPCA.traj = (exTrajStruct(i).avgSmoothFR.traj-mu)*W;
        tbHoldTrajStruct(i).avgdPCA.traj = (tbHoldTrajStruct(i).avgSmoothFR.traj-mu)*W;
        %Get VAF
        %trajStruct(i).avgPCA.VAF = 100.*(diag(cov(allAvgs*coeff))')./totalVar;
        
        %Add all traces to trajStruct; Compute 95% CI
        numExTrials = size(exTrajStruct(i).allSmoothFR,2);
        exAlldPCA = NaN(numPts,size(exTrajStruct(i).avgdPCA.traj,2),numExTrials);
        for j = 1:numExTrials
          numTrialPts = size(exTrajStruct(i).allSmoothFR(j).traj,1);
          if numTrialPts < numPts
            exAlldPCA(1:numTrialPts,:,j) = (exTrajStruct(i).allSmoothFR(j).traj(1:numTrialPts,:)-mu)*W;  
          else   
            exAlldPCA(:,:,j) = (exTrajStruct(i).allSmoothFR(j).traj(1:numPts,:)-mu)*W;
          end
        end
        
        numHoldTrials = size(tbHoldTrajStruct(i).allSmoothFR,2);
        holdAlldPCA = NaN(numHoldPts,size(tbHoldTrajStruct(i).avgdPCA.traj,2),numHoldTrials);
        for j = 1:numHoldTrials
          numTrialPts = size(tbHoldTrajStruct(i).allSmoothFR(j).traj,1);
          if numTrialPts < numHoldPts
            holdAlldPCA(1:numTrialPts,:,j) = (tbHoldTrajStruct(i).allSmoothFR(j).traj(1:numTrialPts,:)-mu)*W;  
          else   
            holdAlldPCA(:,:,j) = (tbHoldTrajStruct(i).allSmoothFR(j).traj(1:numHoldPts,:)-mu)*W;
          end
          holdAlldPCA(:,:,j) = (tbHoldTrajStruct(i).allSmoothFR(j).traj(1:numHoldPts,:)-mu)*W;
        end

        %Get confidence intervals
        exTrajStruct(i).avgdPCA.CI = std(exAlldPCA,0,3)/sqrt(numExTrials);
        %tbHoldTrajStruct(i).avgdPCA.CI = 1.96.*std(holdAlldPCA,0,3)/sqrt(numHoldTrials);
        tbHoldTrajStruct(i).avgdPCA.CI = std(holdAlldPCA,0,3)/sqrt(numHoldTrials);
    end
    
%% Plot
    % Plot Marginal dims
    postureMargID = find(strcmp(margNames,'P'));
    targetMargID = find(strcmp(margNames,'T'));
    CIMargID = find(strcmp(margNames,'CI'));
    intMargID = find(strcmp(margNames,'PTI'));
    
    postureDims = find(whichMarg==postureMargID,2);
    targetDims = find(whichMarg==targetMargID,2);
    CIDims = find(whichMarg==CIMargID,2);
    intDims = find(whichMarg==intMargID,2);
    
    dimList = [postureDims(1),targetDims(1),CIDims(1),postureDims(2),targetDims(2),CIDims(2)];

    plotTargetList = [1,5];
        
    timeOffset = 150;
    f = figure; fs = 14;
    f.Position = [200,200,1225,464];
    postureInd = 1;
    for posture = postureList
        for target = plotTargetList
            holdTraj = tbHoldTrajStruct([tbHoldTrajStruct.posture]==posture & [tbHoldTrajStruct.target]==target).avgdPCA.traj(1:numHoldPts,:);
            holdCI = tbHoldTrajStruct([tbHoldTrajStruct.posture]==posture & [tbHoldTrajStruct.target]==target).avgdPCA.CI(1:numHoldPts,:);
            holdTime = tbHoldTrajStruct([tbHoldTrajStruct.posture]==posture & [tbHoldTrajStruct.target]==target).avgSmoothFR.timestamps(1:numHoldPts);
            
            exTraj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgdPCA.traj(1:numPts,:);
            exCI = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgdPCA.CI(1:numPts,:);
            exTime = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.timestamps(1:numPts);
            exTime = exTime + timeOffset; 
            
            dimInd = 1;
            for dim = dimList
                subplot(2,3,dimInd); hold on
                if target == 1 || target == 3

                    shadedErrorBar(holdTime,holdTraj(:,dim),holdCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                    hold on;
                    shadedErrorBar(exTime,exTraj(:,dim),exCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                    
                elseif target == 5 || target == 7
                    %plot(time,traj(:,dim),'--','Color',pcmap(posture,:),'LineWidth',2);
                    
                    shadedErrorBar(holdTime,holdTraj(:,dim),holdCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)})
                    hold on; 
                    shadedErrorBar(exTime,exTraj(:,dim),exCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                    
                end
                dimInd = dimInd + 1;
            end
            
        end
        postureInd = postureInd + 1;
    end
            
         dimInd = 1;
        for dim = dimList
            subplot(2,3,dimInd); hold on
            ax = gca; xlimits = ax.XLim; ylimits = ax.YLim;
            ax.TickDir = 'out';
            set(gca,'fontname','arial'); set(gca,'fontsize',fs);
            dimVar = explVar.componentVar(dim);
            text(xlimits(1)+10,ylimits(2)*0.9,[num2str(round(dimVar,1)),'%'],'FontSize',fs,'FontName','arial')
            %yticklabels({})
            if ismember(dimInd,[1,2,3])
                xticklabels({});
            end
            if ismember(dimInd,[1,2,3])
                xticks([-100,0,150,250])
                xticklabels({});
            end
            if ismember(dimInd,[4,5,6])
                xticks([-100,0,150,250])
            end 
            if dimInd == 5
                xlabel('time (ms)')
            end 
            if dimInd == 1
                title('Posture');
            end
            if dimInd == 2
                title('Target');
            end
            if dimInd == 3
                title('Time');
            end
            dimInd = dimInd + 1;
        end
        
        sgtitle(dataset);
    
    %Single Stacked Version
    f = figure; fs = 14;
    postureInd = 1;
    for posture = postureList
        for target = plotTargetList
            holdTraj = tbHoldTrajStruct([tbHoldTrajStruct.posture]==posture & [tbHoldTrajStruct.target]==target).avgdPCA.traj(1:numHoldPts,:);
            holdCI = tbHoldTrajStruct([tbHoldTrajStruct.posture]==posture & [tbHoldTrajStruct.target]==target).avgdPCA.CI(1:numHoldPts,:);
            holdTime = tbHoldTrajStruct([tbHoldTrajStruct.posture]==posture & [tbHoldTrajStruct.target]==target).avgSmoothFR.timestamps(1:numHoldPts);
            
            exTraj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgdPCA.traj(1:numPts,:);
            exCI = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgdPCA.CI(1:numPts,:);
            exTime = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.timestamps(1:numPts);
            exTime = exTime + timeOffset; 
            
            dimInd = 1;
            for dim = dimList(1:3)
                subplot(3,1,dimInd); hold on
                if target == 1 || target == 3

                    shadedErrorBar(holdTime,holdTraj(:,dim),holdCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                    hold on;
                    shadedErrorBar(exTime,exTraj(:,dim),exCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                    
                elseif target == 5 || target == 7
                    %plot(time,traj(:,dim),'--','Color',pcmap(posture,:),'LineWidth',2);
                    
                    shadedErrorBar(holdTime,holdTraj(:,dim),holdCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)})
                    hold on; 
                    shadedErrorBar(exTime,exTraj(:,dim),exCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                    
                end
                dimInd = dimInd + 1;
            end
            
        end
        postureInd = postureInd + 1;
    end
    
    
    dimInd = 1;
    for dim = dimList(1:3)
        subplot(3,1,dimInd); hold on
        ax = gca; xlimits = ax.XLim; ylimits = ax.YLim;
        ax.TickDir = 'out';
        set(gca,'fontname','arial'); set(gca,'fontsize',fs);
        dimVar = explVar.componentVar(dim);
        text(xlimits(1)+10,ylimits(2)*0.9,[num2str(round(dimVar,1)),'%'],'FontSize',fs,'FontName','arial')
        %yticklabels({})
        if ismember(dimInd,[1,2])
            xticklabels({});
        end
        if dimInd == 3
            xticks([-100,0,150,250])
            xlabel('time (ms)')
        end 
        if dimInd == 1
            ylabel('Posture (SD)');
        end
        if dimInd == 2
            ylabel('Target (SD)');
        end
        if dimInd == 3
            ylabel('Time (SD)');
        end
        dimInd = dimInd + 1;
    end

    
    
    
    
    
    %Big stacked version 
    f = figure; fs = 14;
    postureInd = 1;
    for posture = postureList
        for target = plotTargetList
            holdTraj = tbHoldTrajStruct([tbHoldTrajStruct.posture]==posture & [tbHoldTrajStruct.target]==target).avgdPCA.traj(1:numHoldPts,:);
            holdCI = tbHoldTrajStruct([tbHoldTrajStruct.posture]==posture & [tbHoldTrajStruct.target]==target).avgdPCA.CI(1:numHoldPts,:);
            holdTime = tbHoldTrajStruct([tbHoldTrajStruct.posture]==posture & [tbHoldTrajStruct.target]==target).avgSmoothFR.timestamps(1:numHoldPts);
            
            exTraj = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgdPCA.traj(1:numPts,:);
            exCI = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgdPCA.CI(1:numPts,:);
            exTime = exTrajStruct([exTrajStruct.posture]==posture & [exTrajStruct.target]==target).avgSmoothFR.timestamps(1:numPts);
            exTime = exTime + timeOffset; 
            
            dimInd = 1;
            for dim = [1,5,2,6,3]
                subplot(5,1,dimInd); hold on
                if target == 1 || target == 3

                    shadedErrorBar(holdTime,holdTraj(:,dim),holdCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                    hold on;
                    shadedErrorBar(exTime,exTraj(:,dim),exCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                    
                elseif target == 5 || target == 7
                    %plot(time,traj(:,dim),'--','Color',pcmap(posture,:),'LineWidth',2);
                    
                    shadedErrorBar(holdTime,holdTraj(:,dim),holdCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)})
                    hold on; 
                    shadedErrorBar(exTime,exTraj(:,dim),exCI(:,dim),'lineprops',{'LineWidth',1.5,'Color',pcmap(posture,:)});
                    
                end
                dimInd = dimInd + 1;
            end
            
        end
        postureInd = postureInd + 1;
    end
    
    
    
    
    
    
% plot_dpca_timecourses(exTrajStruct,numPts,postureList,plotTargetList,W,dimList,pcmap,explVar,dataset)
% 
% plot_dpca_timecourses(tbHoldTrajStruct,numPts,postureList,plotTargetList,W,dimList,pcmap,explVar,dataset)     

%         if saveFig
%             saveas(gcf,fullfile(saveDir,[dataset,'_dpcaTraj.svg']));
%         end

% %% Local plotting function 
% function [] = plot_dpca_timecourses(trajStruct,numPts,postureList,plotTargetList,W,dimList,pcmap,explVar,dataset)
%         f = figure; f.Position = [200 200 850 400]; fs = 14;
%         postureInd = 1;
%         time = trajStruct(1).avgSmoothFR.timestamps(1:numPts);
%         for posture = postureList
%             for target = plotTargetList
%                 traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(1:numPts,:);
%                 traj = traj*W;
%                 dimInd = 1;
%                 for dim = dimList
%                     subplot(2,3,dimInd); hold on
%                     if target == 1 || target == 3
%                         plot(time,traj(:,dim),'Color',pcmap(posture,:),'LineWidth',2);
%                     elseif target == 5 || target == 7
%                         plot(time,traj(:,dim),'--','Color',pcmap(posture,:),'LineWidth',2);
%                     end
%                     dimInd = dimInd + 1;
%                 end
%             end
%             postureInd = postureInd + 1;
%         end
% 
%         dimInd = 1;
%         for dim = dimList
%             subplot(2,3,dimInd); hold on
%             ax = gca; xlimits = ax.XLim; ylimits = ax.YLim;
%             set(gca,'fontname','arial'); set(gca,'fontsize',fs);
%             dimVar = explVar.componentVar(dim);
%             text(xlimits(1)+10,ylimits(2)*0.9,[num2str(round(dimVar,1)),'%'],'FontSize',fs,'FontName','arial')
%             %yticklabels({})
%             if ismember(dimInd,[1,2,3])
%                 xticklabels({});
%             end
%             if dimInd == 5
%                 xlabel('time (ms)')
%             end 
%             if dimInd == 1
%                 title('Posture');
%             end
%             if dimInd == 2
%                 title('Target');
%             end
%             if dimInd == 3
%                 title('Time');
%             end
%             dimInd = dimInd + 1;
%         end
%         
%         sgtitle(dataset);
% end