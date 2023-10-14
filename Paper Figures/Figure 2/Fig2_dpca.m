clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 2';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set parameters
    bciAnalysisWindow = [];
    isoAnalysisWindow = [];
    reachAnalysisWindow = [];

%% Main Loop    
    %{'E20200317','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}
    for datasetList = {'E20200316'}%{'E20200317','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}%{'R20200221'}%{'E20200317','E20200116','E20210706'}       
        %% Load Data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        %% Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318'}
                %trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            case {'N20171215','N20180221'}
                %trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
            case {'R20201020','R20201021'}
                %trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
            case 'E20210901'
                taskID = [Data.conditionData]; taskID = [taskID.taskID];
                Data = Data(taskID==1);
                %trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
            %Iso
            case 'E20200116'
                %trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',500}};
            %Reaching
            case 'E20210706'
                %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
            case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
                %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
            case{'R20200221','R20200222'}
                %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
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
        switch dataset
            case {'E20200316'}
                pcmap = flip(pcmap);
            case {'N20171215'}     
                pcmap = vertcat(pcmap(3,:),pcmap(1,:),pcmap(5,:));
            case {'E20210901'}
                pcmap = vertcat(pcmap(1,:),pcmap(3,:),pcmap(3,:),pcmap(5,:),pcmap(4,:));
                pcmap = vertcat(pcmap(4,:),pcmap(5,:),pcmap(5,:),pcmap(1,:),pcmap(3,:));
            case {'R20201020'}
                pcmap = vertcat(pcmap(2,:),pcmap(1,:));
        end
        
        %Form X
        [X] = getX(trajStruct,minNumTimestamps,postureList,numPostures,targetList,numTargets,numChannels);   
        
        %Do dPCA
        Xdpca = permute(X,[4,3,2,1]);

        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        margNames = {'P', 'T', 'CI', 'PTI'};
        margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

        N = numChannels; P = numPostures; T = numTargets;
        timePts = minNumTimestamps;
        time = trajStruct(1).avgZSmoothFR.timestamps(1:minNumTimestamps);
        timeEvents = [];

        [W,V,whichMarg] = dpca(Xdpca, 10, ...
            'combinedParams', combinedParams);

        explVar = dpca_explainedVariance(Xdpca, W, V, ...
            'combinedParams', combinedParams);

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
            %Get VAF
            trajStruct(i).PTOrth.VAF =  [explVar.componentVar(postureDims(1)),explVar.componentVar(targetDims(1)),explVar.componentVar(targetDims(2))];
            %100.*(diag(cov(allTraj*PTOrth))')./totalVar;
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
        %clf; close all
    end