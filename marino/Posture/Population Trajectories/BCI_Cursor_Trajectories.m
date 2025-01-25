clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\Figure 1';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
    %% Run loop for each dataset    
    for datasetList = {'E20200317'} 
        %datasetList = {'E20200317','N20180221'} 
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'bciCursorTraj'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            case 'E20200317'
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            case 'N20180221'
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
        end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numConditions = size(trajStruct,2);
    
        %Get pcmap for number of postures
        pcmap = orli(1:numPostures,:);
        lightpcmap = rgb2hsv(pcmap);
        lightpcmap(:,2)=.2;
        lightpcmap =hsv2rgb(lightpcmap);
        
        %Create plots
        fs = 14;
        for posture = postureList
            figure; hold on
            for target = targetList
                allTraj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allBciCursorTraj;
                %Get 10 trajectories per target 
                numTraj = size(allTraj,2);
                trajIDs = randsample(numTraj,10);
                for trajID = trajIDs'
                    traj = allTraj(trajID).traj;
                    plot(traj(:,1),traj(:,2),'Color',lightpcmap(posture,:))
                end
            end
            for target = targetList
               avgTraj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgBciCursorTraj.traj; 
               plot(avgTraj(:,1),avgTraj(:,2),'LineWidth',3,'Color',pcmap(posture,:)) 
            end
            xlim([-125 125])
            ylim([-125 125])
            xticks([-100 0 100])
            yticks([-100 0 100])    
            xlabel('x (mm)')
            ylabel('y (mm)')
            set(gca,'fontname','arial')
            set(gca,'fontsize',fs)
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'p',num2str(posture),'_BCI_cursorTraj.svg']));
            end
        end

    end