clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\Consistent Posture Signal';
    set(0, 'DefaultFigureRenderer', 'painters');
    fs = 14;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;
    
%% Load Data, subselect
    dataset = 'E20200314';
    [Data,zScoreParams] = loadData(dataset);
    
    %Get trajStruct
    binWidth = 25; kernelStdDev = 25;
    trialInclStates = struct('trialName','','inclStates',[]);
    trajFields = {'zSmoothFR','markerVel','marker'};
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    
%% Get trajStruct   
    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
        trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',25},{'state','Step 2','first',0}};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        %trialInclStates(2).inclStates = {{'state','Reach','first',25},{'state','Reach','first',300}};
        trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-500},{'state','Hold','first',0}};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
   
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

    taskIDs = struct('ID',[],'task','');
    taskIDs(1).ID = 1; taskIDs(1).task = 'BC';
    taskIDs(2).ID = 2; taskIDs(2).task = 'HC';
    taskIDs(3).ID = 3; taskIDs(3).task = 'Iso';
    taskIDs(4).ID = 4; taskIDs(4).task = 'All';
   
%% Examine reach kinematics on individual trials
    trial = 574;
    trialData = Data(trial);
    trialInclStates = struct('trialName','','inclStates',[]);
    trialInclStates(1).trialName = {'HC_CenterOut_ForceBar_20200314'};
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-500},{'state','Hold','first',0}};
    %trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Hold','first',0}};
    [markerVel,velTime] = getStatesTraj20220419(trialData,trialInclStates,'markerVel',binWidth,kernelStdDev,'timeRelToTrialStart',true);
    [markerPos,posTime] = getStatesTraj20220419(trialData,trialInclStates,'marker',binWidth,kernelStdDev,'timeRelToTrialStart',true);

    figure
    plot(velTime,markerVel(:,1))
    hold on
    plot(velTime,markerVel(:,2))


    figure
    plot(markerPos(:,1),markerPos(:,2))
    hold on
    targetLoc = Data(trial).targetData.targetLoc;
    targetSize = Data(trial).targetData.targetSize;
    circle(targetLoc(1,1),targetLoc(1,2),targetSize,'k')
    axis equal

%% Examine force kinematics on individual trials
 trial = 318;
    trialData = Data(trial);
    trialInclStates = struct('trialName','','inclStates',[]);
      trialInclStates(1).trialName = {'IsometricForce_1D'};
        trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',300}};
   
    [force,forceTime] = getStatesTraj20220419(trialData,trialInclStates,'force',binWidth,kernelStdDev,'timeRelToTrialStart',false);
    [forceCursor,forceCursorTime] = getStatesTraj20220419(trialData,trialInclStates,'forceCursor',binWidth,kernelStdDev,'timeRelToTrialStart',false);

    figure; 
    subplot(4,1,1); plot(forceTime,force(:,1))
    hold on
    subplot(4,1,2); plot(forceTime,force(:,2))
    subplot(4,1,3); plot(forceTime,force(:,3))
    subplot(4,1,4); plot(forceTime,force(:,4))

    figure
    plot(forceCursorTime,forceCursor(:,2))
    hold on
    targetLoc = trialData.targetData.targetLoc(1,2);
    ax = gca;
    xlim = ax.XLim;
    plot([xlim(1) xlim(2)],[targetLoc targetLoc],'r-')
    %targetSize = Data(trial).targetData.targetSize;
    %circle(targetLoc(1,1),targetLoc(1,2),targetSize,'k')


%% Plot traces of each trial type; select appropriate time ranges of neural activity
 %Plot Neuron Time Course
    task = 3;
    tempTrajStruct = trajStruct([trajStruct.task]==task);
    targetList = unique([tempTrajStruct.target]);
    f = figure;
    f.Position = [20 20 900 500];
    for posture = 1
        for target = targetList
            time = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).avgSmoothFR.timestamps;
            traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).avgSmoothFR.traj;
            numNeurons = size(traj,2);
            for i = 1:25
               subplot(5,5,i)
                   plot(time,traj(:,i),'Color',tcmap(target,:));
                   hold on
                   xlabel('t (ms)')
                   ylabel('FR (std)')
            end
        end
    end
%     sgtitle(epoch{1,1})
%     if saveFig
%         saveas(gcf,fullfile(saveDir,[dataset,'_',epoch{1,1},'_NeuronTimeCourse.svg']));
%     end