clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'D:\Posture Paper\20221123\Figures\Supplement\All Task dPCA';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Run Loop    
    for datasetList = {'E20210901'}%{'E20200317','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}%{'R20200221'}%{'E20200317','E20200116','E20210706'}
       clf; close all
        %{'E20200317','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}
        % Load Data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318'}
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
            case 'E20210901'
                taskID = [Data.conditionData]; taskID = [taskID.taskID];
                Data = Data(taskID==1);
                trialInclStates(1).trialName = {'BCI Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
            %Iso
            case 'E20200116'
                trialInclStates(1).trialName = {'IsometricForce_1D'};   
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',250}};
            %Reaching
            case 'E20210706'
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
            case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
            case{'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
        end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        
        %Keep postures 1-4 only for earl
        switch dataset
            case {'E20210706'}
                    trajStruct = trajStruct(ismember([trajStruct.posture],[1,2,3,4]));
        end
        

        % Get minimum number of timestamps in condition averages
        numTimestamps = [];
        for i = 1:size(trajStruct,2)
            numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
        end
        [minNumTimestamps,i] = min(numTimestamps);

        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);

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

        Xdpca = permute(X,[4,3,2,1]);

        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        margNames = {'P', 'T', 'CI', 'PTI'};
        margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

        N = numChannels;
        P = numPostures;
        T = numTargets;
        timePts = minNumTimestamps;
        time = trajStruct(1).avgSmoothFR.timestamps(1:minNumTimestamps);
        timeEvents = [];

        [W,V,whichMarg] = dpca(Xdpca, 20, ...
            'combinedParams', combinedParams);

        explVar = dpca_explainedVariance(Xdpca, W, V, ...
            'combinedParams', combinedParams);

        dpca_plot(Xdpca, W, V, @dpca_plot_default, ...
            'explainedVar', explVar, ...
            'marginalizationNames', margNames, ...
            'marginalizationColours', margColours, ...
            'whichMarg', whichMarg,                 ...
            'time', time,                        ...
            'timeEvents', timeEvents,               ...
            'timeMarginalization', 3, ...
            'legendSubplot', 1);

        % Plot Variance Pie
        fs = 14;
        figure; set(gca,'fontname','arial'); set(gca,'fontsize',fs);
        
        p = pie(explVar.totalMarginalizedVar/explVar.totalVar);
        txt = {'Posture: ';'Target: ';'CI: ';'Interaction: '}; 

        pText = findobj(p,'Type','text');
        percentValues = get(pText,'String'); 

        combinedtxt = strcat(txt,percentValues); 

        for i = 1:size(pText,1)
            pText(i).String = combinedtxt(i);
            pText(i).FontSize = fs;
            pText(i).FontName = 'arial';
        end

        title(dataset);
        
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_dpcaPieChart.svg']));
        end

        % Plot Marginal dims
        postureMargID = find(strcmp(margNames,'P'));
        targetMargID = find(strcmp(margNames,'T'));
        CIMargID = find(strcmp(margNames,'CI'));

        postureDims = find(whichMarg==postureMargID,2);
        targetDims = find(whichMarg==targetMargID,2);
        CIDims = find(whichMarg==CIMargID,2);
        dimList = [postureDims(1),targetDims(1),CIDims(1),postureDims(2),targetDims(2),CIDims(2)];

        plotTargetList = [1,5];
        if strcmpi(dataset,'E20200116')
            plotTargetList = [3,7];
        end
        
        f = figure; f.Position = [200 200 850 400];
        postureInd = 1;
        time = trajStruct(1).avgSmoothFR.timestamps(1:minNumTimestamps);
        for posture = postureList
            for target = plotTargetList
                traj = Xdpca(:,find(postureList==posture),find(targetList==target),:);
                traj = squeeze(traj);
                traj = traj'*W;
                dimInd = 1;
                for dim = dimList
                    subplot(2,3,dimInd); hold on
                    if target == 1 || target == 3
                        plot(time,traj(:,dim),'Color',pcmap(posture,:),'LineWidth',2);
                    elseif target == 5 || target == 7
                        plot(time,traj(:,dim),'--','Color',pcmap(posture,:),'LineWidth',2);
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
            set(gca,'fontname','arial'); set(gca,'fontsize',fs);
            dimVar = explVar.componentVar(dim);
            text(xlimits(1)+10,ylimits(2)*0.9,[num2str(round(dimVar,1)),'%'],'FontSize',fs,'FontName','arial')
            yticklabels({})
            if ismember(dimInd,[1,2,3])
                xticklabels({});
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
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_dpcaTraj.svg']));
        end
    end