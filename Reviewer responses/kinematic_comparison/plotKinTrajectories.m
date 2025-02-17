clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 3';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(5); 
    
%% Load data, getTrajStruct
    dataset = 'R20201021';
    %Load data
    [Data,zScoreParams] = loadData(dataset);
    [Data] = removeShortBCIandIsoTrials(Data,dataset);
    %Get trajStruct
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParamsKinematicComparison(dataset);
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',false);  
    %Remove any conditions for which there weren't enough trials
    [numCondTrials] = getNumCondTrials(trajStruct,'showHist',false);            
    %get trajStructDims
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
    %Get minimum number of condition trials and timestamps
    [minNumCondTrials] = getMinNumCondTrials(trajStruct);
    [minNumTimestamps] = getMinNumTimestamps(trajStruct); 

    
%% Plot Positions using trajStruct 
    for posture = postureList
        figure; hold on
        for target = targetList
            tempTrajStruct = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target);
            if ~isempty(tempTrajStruct)

                % Plot individual trial trajectories (faded)
                numTrials = size(tempTrajStruct.allBciCursorTraj,2);
                for trial = 1:numTrials
                    traj = tempTrajStruct.allBciCursorTraj(trial).traj;
                    %plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
                    plot(traj(:,1), traj(:,2), ...
                     'Color', [tcmap(target,:) 0.2], ... 
                     'LineWidth', 0.5);
                end
                
                % Plot average trajectory (bold)                
                traj = tempTrajStruct.avgBciCursorTraj.traj;                
                plot(traj(:,1),traj(:,2),'Color',tcmap(target,:), 'LineWidth', 2);
                
                
            end
        end
        title(['Posture ',num2str(posture)])
    end