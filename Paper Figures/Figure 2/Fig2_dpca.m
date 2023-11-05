clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 2';
    set(0, 'DefaultFigureRenderer', 'painters');
    
% To make nicest possible projection, use E20200316; go from Step 1 to Step
% 2; then choose first 8 timestamps from each condition avg. No
% regularization. Could try removing really short trials to aid in
% regularization? CURRENTLY YOU'RE GETTING PART OF STEP 2 FOR REALLY SHORT
% TRIALS

%% Set parameters
%     bciAnalysisWindow = [25,225];
%     isoAnalysisWindow = [];
%     reachAnalysisWindow = [-200,300];
     numDPCs = 15;
%     
%% Main Loop  
%{'R20200221'}%{'E20200317','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}%{'R20200221'}%{'E20200317','E20200116','E20210706'}       
    %{'E20200318','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}
    for datasetList = {'E20200318'}%{'E20200318','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}
        %% Load Data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        
%         %% Remove trials less that were shorter than duration of analysis window
%         kinData = [Data.kinData];
%         moveTime = [kinData.moveTime];
%         Data = Data(moveTime >= bciAnalysisWindow(2));
%         figure
%         histogram(moveTime);
        
        
        %% Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        switch dataset
%             %BCI
             case {'E20200316','E20200317','E20200318'}
                %trialInclStates(1).inclStates = {{'state','Step 1','first',50},{'state','Step 1','first',225}};
%             case {'N20171215','N20180221'}
%                 %trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
%             case {'R20201020','R20201021'}
%                 %trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
%             case 'E20210901'
%                 taskID = [Data.conditionData]; taskID = [taskID.taskID];
%                 Data = Data(taskID==1);
%                 %trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
%             %Iso
            case 'E20200116'
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',250}};
%             %Reaching
%             case 'E20210706'
%                 trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',reachAnalysisWindow(1)},{'kin','moveOnsetTime','first',reachAnalysisWindow(2)}};
%                 %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
%             case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
%                 trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',reachAnalysisWindow(1)},{'kin','moveOnsetTime','first',reachAnalysisWindow(2)}};
%                 %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
%             case{'R20200221','R20200222'}
%                 %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',reachAnalysisWindow(1)},{'kin','moveOnsetTime','first',reachAnalysisWindow(2)}};
%                 %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
        end     
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);
    
        %% Keep postures 1-4 only for earl reaching
        if strcmpi(dataset,'E20210706')
            trajStruct = trajStruct(ismember([trajStruct.posture],[1,2,3,4]));
        end
              
        %% Get traj struct dimensions
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 
        
        %% Setup colormap (based on number of postures)
        [pcmap,tcmap,rainbow] = getColorMaps(numPostures);  
%         switch dataset
%             case {'E20200316'}
%                 pcmap = flip(pcmap);
%             case {'N20171215'}     
%                 pcmap = vertcat(pcmap(3,:),pcmap(1,:),pcmap(5,:));
%             case {'E20210901'}
%                 pcmap = vertcat(pcmap(1,:),pcmap(3,:),pcmap(3,:),pcmap(5,:),pcmap(4,:));
%                 pcmap = vertcat(pcmap(4,:),pcmap(5,:),pcmap(5,:),pcmap(1,:),pcmap(3,:));
%             case {'R20201020'}
%                 pcmap = vertcat(pcmap(2,:),pcmap(1,:));
%         end
%         
        %% Form X and Xfull
        [X] = getX(trajStruct,'avgZSmoothFR',minNumTimestamps,postureList,numPostures,targetList,numTargets);    
        [Xfull] = getXFull(trajStruct,'allZSmoothFR',minNumTimestamps,postureList,numPostures,targetList,numTargets);    
        
        %Permute for use w dPCA (to follow along w examples)
        Xdpca = permute(X,[4,3,2,1]);
        XdpcaFull = permute(Xfull,[4,3,2,1,5]);
        
        %% dPCA Parameters
        
        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        margNames = {'P', 'T', 'CI', 'PTI'};
        margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
        N = numChannels; P = numPostures; T = numTargets;
        time = trajStruct(1).avgZSmoothFR.timestamps(1:minNumTimestamps);
        timeEvents = [];

        %% Regularization 
        %Set up 'numOfTrials' for use with dpca_optimizeLambda
        numOfTrials = zeros(size(Xdpca,1),size(Xdpca,2),size(Xdpca,3));
        for i = 1:size(Xdpca,1)
            for j = 1:size(Xdpca,2)
                for k = 1:size(Xdpca,3)
                    numOfTrials(i,j,k) = sum(~isnan(XdpcaFull(i,j,k,1,:)));
                end
            end
        end
        
        %Compute optimal lambda
        optimalLambda = dpca_optimizeLambda(Xdpca, XdpcaFull, numOfTrials, ...
            'combinedParams', combinedParams, ...
            'simultaneous', true, ...
            'numRep', 10, ...  % increase this number to ~10 for better accuracy
            'filename', 'tmp_optimalLambdas.mat');

        %Compute Cnoise
        Cnoise = dpca_getNoiseCovariance(Xdpca, ...
            XdpcaFull, numOfTrials, 'simultaneous', true);

        %% Compute dPCs
%         [W,V,whichMarg] = dpca(Xdpca, numDPCs, ...
%             'combinedParams', combinedParams, ...
%             'lambda', optimalLambda, ...
%             'Cnoise', Cnoise);

        [W,V,whichMarg] = dpca(Xdpca, numDPCs, ...
            'combinedParams', combinedParams);
        
        explVar = dpca_explainedVariance(Xdpca, W, V, ...
            'combinedParams', combinedParams);
        cumulativeDPCA = explVar.cumulativeDPCA;
        additionalVarExpl = [cumulativeDPCA(1) diff(cumulativeDPCA)];
        
%         dpca_plot(Xdpca, W, V, @dpca_plot_default, ...
%             'explainedVar', explVar, ...
%             'marginalizationNames', margNames, ...
%             'marginalizationColours', margColours, ...
%             'whichMarg', whichMarg,                 ...
%             'time', time,                        ...
%             'timeEvents', timeEvents,               ...
%             'timeMarginalization', 3, ...
%             'legendSubplot', 1);

        % Plot Marginal dims
        postureMargID = find(strcmp(margNames,'P'));
        targetMargID = find(strcmp(margNames,'T'));
        CIMargID = find(strcmp(margNames,'CI'));

        postureDims = find(whichMarg==postureMargID,2);
        targetDims = find(whichMarg==targetMargID,2);
        CIDims = find(whichMarg==CIMargID,2);
        %dimList = [postureDims(1),targetDims(1),CIDims(1),postureDims(2),targetDims(2),CIDims(2)];

        %Orthonormalize pdpca1 and tdpca1-2
        pDPCA = W(:,postureDims);
        tDPCA = W(:,targetDims);
        [PTOrth,~] = qr([pDPCA(:,1),tDPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
             
        
        %Add Projections to trajStruct
        allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
        totalVar = trace(cov(allTraj));
        for i = 1:size(trajStruct,2)
            %Add average traces to trajStruct
            trajStruct(i).PTOrth.traj = trajStruct(i).avgZSmoothFR.traj*PTOrth;
            trajStruct(i).dPCA.traj = trajStruct(i).avgZSmoothFR.traj*W;
            %Get VAF
            trajStruct(i).PTOrth.VAF =  [additionalVarExpl(postureDims(1)),additionalVarExpl(targetDims(1)),additionalVarExpl(targetDims(2))];
            trajStruct(i).dPCA.VAF = explVar.componentVar;
        end

        %Plot - Orthographic
        fs = 14;
        figure
        hold on
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1;
        for posture = postureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                    plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                end
            end
        end
             
        
        grid on
        xticklabels({}); yticklabels({}); zticklabels({}); 
        VAF = round(trajStruct(1).PTOrth.VAF);

        xlabel(['Goal Dim 1 (',num2str(VAF(xDim)),'%)'])
        ylabel(['Goal Dim 2 (',num2str(VAF(yDim)),'%)']) 
        zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
        view([20 10])
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_PTOrth.fig']));
        end
        
        
        %Plot vs time
        figure 
        timePts = 1:minNumTimestamps;
        for posture = postureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).dPCA.traj(timePts,:); 
                    time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgZSmoothFR.timestamps(timePts); 
                    for dim = 1:10
                       subplot(2,5,dim)
                       plot(time,traj(:,dim),'Color',pcmap(posture,:),'LineWidth',2);
                       hold on;                        
                    end
                end
            end
        end
          
        %Format
        maxYRange = 0;
        for dim = 1:10
            subplot(2,5,dim)
            ax = gca;
            ylimits = ax.YLim;
            yRange = ylimits(2)-ylimits(1);
            if yRange > maxYRange
                maxYRange = yRange;
            end
        end
        for dim = 1:10
            VAF = round(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).dPCA.VAF);
            subplot(2,5,dim);
            ax = gca;
            ylimits = ax.YLim;
            yMid = ylimits(1) + 0.5*(ylimits(2)-ylimits(1));
            ylim([yMid-maxYRange/2 yMid+maxYRange/2])
            ylabel(['dPC ',num2str(dim),' (',num2str(VAF(dim)),'%)'])
            yticks([min(ax.YTick) max(ax.YTick)]);
            ax.TickDir = 'out';
            if dim < 4
                xticks([])
            else
                xticks([-200 0 300]);
            end
            if dim == 5
                %xlabel('time (ms)')
            end
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
        end           

        
        clf; close all
    end