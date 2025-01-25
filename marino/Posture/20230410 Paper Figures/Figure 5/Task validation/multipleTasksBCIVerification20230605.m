clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\Analysis window verification\BCI Validation';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Load data; get trajStruct
    multipleTasksDatasetList = {'E20200311','E20200312','E20200313','E20200314'};
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};

    for datasetList = {'R20201020','R20201021'}
        %Load data
         dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        %Set up common part of trialInclStates       
        binWidth = 25; kernelStdDev = 25;
        trialInclStates = struct('trialName','','inclStates',[]);
        trajFields = {'zSmoothFR'};    
                
        %Keep only BCI 
        switch dataset
            case {'E20200311','E20200312','E20200313','E20200314'}
                conditionData = [Data.conditionData]; taskID = [conditionData.taskID];
                Data = Data([taskID==1]);
                condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
                trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',-100},{'state','Step 2','first',0}};
            %BCI
            case {'E20200316','E20200317','E20200318','E20200319'}
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',-100},{'state','Step 2','first',0}};
                task = 'bci';
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
                task = 'bci';
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
                task = 'bci';
        end
        %Get trajStruct
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        %Get kinStruct
        kinFields = {'moveTime'};
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
        case {'E20200316','E20200317','E20200318','E20200319'}
            postureList = [1,5];
        case {'N20171215','N20180221','R20201020','R20201021'}
            postureList = [1,2];
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

        
        %movement time histogram
        figure;
        plotMat = [];
        plotMatInd = 1;
        grp = [];
        grpInd = 1;
        for target = targetList
            for posture = postureList
                allMoveTime = kinStruct([kinStruct.target]==target & [kinStruct.posture]==posture).allMoveTime;
                plotMat(plotMatInd:plotMatInd+length(allMoveTime)-1,1) = allMoveTime;
                plotMatInd = plotMatInd + length(allMoveTime);
                grp = [grp,(grpInd-1)*ones(1,length(allMoveTime))];
                grpInd = grpInd + 1;
            end
        end
        boxplot(plotMat,grp); 
        ylabel('move time (ms)')
        
        xticklabels({...
            ['P',num2str(postureList(1)),'T',num2str(targetList(1))],...
            ['P',num2str(postureList(2)),'T',num2str(targetList(1))],...
            ['P',num2str(postureList(1)),'T',num2str(targetList(2))],...
            ['P',num2str(postureList(2)),'T',num2str(targetList(2))]});
        title(num2str(dataset))
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_MoveTime.svg']));
        end   
        
        close all;
    end