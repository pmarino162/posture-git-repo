clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\Analysis window verification\Reaching Validation';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Load data; get trajStruct
    multipleTasksDatasetList = {'E20200311','E20200312','E20200313','E20200314'};
    reachingDatasetList = {'E20210706','E20210707','E20210708','E20210709','E20210710',...
        'N20190222','N20190226','N20190227','N20190228','N20190307',...
        'R20200221','R20200222'};
    for datasetList = reachingDatasetList
        %Load data
         dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        %Set up common part of trialInclStates       
        binWidth = 25; kernelStdDev = 25;
        trialInclStates = struct('trialName','','inclStates',[]);
        trajFields = {'zSmoothFR','marker','markerVel'};    
        
        
        %Keep only iso force
        switch dataset
            case {'E20200311','E20200312','E20200313','E20200314'}
                conditionData = [Data.conditionData]; taskID = [conditionData.taskID];
                Data = Data([taskID==2]);
                condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).trialName = {'HC_CenterOut_ForceBar_20200314'};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-275},{'kin','moveOnsetTime','first',250}};
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).trialName = {'GridReaching'};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-275},{'kin','moveOnsetTime','first',250}};
            case {'N20190222','N20190226','N20190227','N20190228','N20190307'}
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',250}};
            case {'R20200221','R20200222'}
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',250}};
        end
        %Get trajStruct
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        %Get kinStruct
        kinFields = {'rxnTime','reachTime'};
        [kinStruct] = getKinStruct(Data,condFields,kinFields);
        
    %% Do PCA on neural activity; add to trajStruct
        numManPCs = 10;
        %Get minimum number of trials and timestamps
        numTimestamps = [];
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allSmoothFR,2);
           numTimestamps = [numTimestamps,size(trajStruct(i).avgSmoothFR.timestamps,2)]; 
        end
        [minNumTimestamps,i] = min(numTimestamps);
        numPts = minNumTimestamps;

        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);

        allTraj = NaN(numConditions*numPts,numChannels);
        j = 1;
        for i = 1:numConditions
           allTraj(j:j+numPts-1,:) = trajStruct(i).avgSmoothFR.traj(1:numPts,:);
           j = j + numPts;
        end

        [allPCs,~,~,~,explained,allMu] = pca(allTraj);
        for i = 1:size(trajStruct,2)
           trajStruct(i).avgPCTraj.traj = (trajStruct(i).avgSmoothFR.traj-allMu)*allPCs(:,1:numManPCs); 
        end



    %% Plot mean force and neural activity for two posture and targets
    targetList = [1,5];
    switch dataset
        case {'E20200311','E20200312','E20200313','E20200314'}
            postureList = [1,3];
        case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
            postureList = [2,6];
            pcmap = vertcat(pcmap,pcmap);
        case {'N20190222','N20190226','N20190227','N20190228','N20190307','R20200221','R20200222'}
            postureList = [1,2];
    end
    
        %Marker Pos and Vel
        figure;
        for target = targetList
            for posture = postureList
                pos = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj;
                posTime = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.timestamps;
                vel = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarkerVel.traj;
                velTime = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarkerVel.timestamps;
                for dim = 1:4
                    subplot(4,1,dim)
                    switch dim
                        case 1
                            traj = pos(:,1);
                            time = posTime;
                        case 2
                            traj = pos(:,2);
                            time = posTime;
                        case 3
                            traj = vel(:,1);
                            time = velTime;
                        case 4
                            traj = vel(:,2);
                            time = velTime;
                    end
                    if target == 1 || target == 3
                        plot(time,traj,'LineWidth',1.5,'Color',pcmap(posture,:));
                        hold on;
                    elseif target == 5 || target == 7
                        plot(time,traj,'--','LineWidth',1.5,'Color',pcmap(posture,:));
                        hold on;
                    end
                end
            end
        end
        for dim = 1:4
            subplot(4,1,dim)
            switch dim 
                case 1
                    ylabel('X (mm)')
                case 2
                    ylabel('Y (mm)')
                case 3
                    ylabel('V_X (m/s)')
                case 4
                    ylabel('V_Y (m/s)')
            end
        end
        xlabel('time (ms)')
        sgtitle(num2str(dataset))
        if saveFig
             saveas(gcf,fullfile(saveDir,[dataset,'_HandKin.svg']));
        end

        %PCs
        f = figure; f.Position = [100 100 750 700];
        for target = targetList
            for posture = postureList
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPCTraj.traj;
                time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
                for dim = 1:10
                subplot(5,2,dim)
                    if target == 1 || target == 3
                        plot(time,traj(:,dim),'LineWidth',1.5,'Color',pcmap(posture,:));
                        hold on;
                    elseif target == 5 || target == 7
                        plot(time,traj(:,dim),'--','LineWidth',1.5,'Color',pcmap(posture,:));
                        hold on;
                    end
                end
            end
        end
        for dim = 1:10
            subplot(5,2,dim)
            ylabel(['PC ',num2str(dim),' (',num2str(round(explained(dim))),'%)'])
        end
        xlabel('time (ms)')
        sgtitle(num2str(dataset))
        if saveFig
             saveas(gcf,fullfile(saveDir,[dataset,'_PCs.svg']));
        end

        %Neurons
        numGroups = ceil(numChannels/49);
        for group = 1:numGroups
            channelList = (group-1)*49 + 1:group*49;
            channelList(channelList > numChannels) = [];
            f = figure; f.Position = [0 0 800 800];
            for target = targetList
                for posture = postureList
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(:,channelList);
                    time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps;
                    for dim = 1:size(channelList,2)
                    subplot(7,7,dim)
                        if target == 1 || target == 3
                            plot(time,traj(:,dim),'LineWidth',1.5,'Color',pcmap(posture,:));
                            hold on;
                        elseif target == 5 || target == 7
                            plot(time,traj(:,dim),'--','LineWidth',1.5,'Color',pcmap(posture,:));
                            hold on;
                        end
                    end
                end
            end
            sgtitle(num2str(dataset))
            if saveFig
                 saveas(gcf,fullfile(saveDir,[dataset,'_Neurons_Group_',num2str(group),'.svg']));
            end      
        end

        
        %reaction time and movement time histogram
        figure;
        plotMat = [];
        plotMatInd = 1;
        grp = [];
        grpInd = 1;
        for target = targetList
            for posture = postureList
                allRxnTime = kinStruct([kinStruct.target]==target & [kinStruct.posture]==posture).allRxnTime;
                plotMat(plotMatInd:plotMatInd+length(allRxnTime)-1,1) = allRxnTime;
                plotMatInd = plotMatInd + length(allRxnTime);
                grp = [grp,(grpInd-1)*ones(1,length(allRxnTime))];
                grpInd = grpInd + 1;
            end
        end
        boxplot(plotMat,grp); 
        ylabel('reaction time (ms)')
        xticklabels({...
            ['P',num2str(postureList(1)),'T',num2str(targetList(1))],...
            ['P',num2str(postureList(2)),'T',num2str(targetList(1))],...
            ['P',num2str(postureList(1)),'T',num2str(targetList(2))],...
            ['P',num2str(postureList(2)),'T',num2str(targetList(2))]});
        %xticklabels({'P1T3','P5T3','P1T7','P5T7'})
        title(num2str(dataset))
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_RxnTime.svg']));
        end   
        
        figure;
        plotMat = [];
        plotMatInd = 1;
        grp = [];
        grpInd = 1;
        for target = targetList
            for posture = postureList
                allReachTime = kinStruct([kinStruct.target]==target & [kinStruct.posture]==posture).allReachTime;
                plotMat(plotMatInd:plotMatInd+length(allReachTime)-1,1) = allReachTime;
                plotMatInd = plotMatInd + length(allReachTime);
                grp = [grp,(grpInd-1)*ones(1,length(allReachTime))];
                grpInd = grpInd + 1;
            end
        end
        boxplot(plotMat,grp); 
        ylabel('reach time (ms)')
        
        xticklabels({...
            ['P',num2str(postureList(1)),'T',num2str(targetList(1))],...
            ['P',num2str(postureList(2)),'T',num2str(targetList(1))],...
            ['P',num2str(postureList(1)),'T',num2str(targetList(2))],...
            ['P',num2str(postureList(2)),'T',num2str(targetList(2))]});
        title(num2str(dataset))
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_ReachTime.svg']));
        end   
        
        close all;
    end