clear; clc; clf; close all;

%% Setup Save
    day = '15';
    month = '07';
    year = '2021';
    dateStr = [year,month,day];
    savePath = fullfile('D:\Animals\Earl',year,month,dateStr,'Daily Plots\');
    saveFig = false;
    
%% Load Data
load('D:\Animals\Earl\2021\07\20210715\03_brainControl\Earl20210715_03_brainControl_SI_translated.mat');
%     load('D:\Animals\Earl\2021\07\20210714\03_brainControl\Earl20210714_03_brainControl_SI_translated.mat');
%     load('D:\Animals\Earl\2021\07\20210713\03_brainControl\Earl20210713_03_brainControl_SI_translated.mat');
    Data = getDataStruct(Data,'getForce',true,'forceSetup','shoulder_posture_bci','getKin',true);

%% Get only successful data 
    sucData = Data([Data.trialStatus]==1);

%% Get Reach Trajectories
    condFields = {{'target','targetData','targetID'}};
    trajFields = {'decoderCursorTraj','force'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'CenterOut_ForceBar_BC'};
        trialInclStates(1).inclStates = {'BC Freeze','Touch Bar BC'}; trialInclStates(1).inclOccurrence = {'first','first'};
%         trialInclStates(1).addTimeToBeginning = {0}; trialInclStates(1).addTimeToEnd = {0};
    trajStruct = getTrajStruct(sucData,condFields,trajFields,trialInclStates);    

%% Get Kinematic Data
    condFields = {{'target','targetData','targetID'}};
    kinFields = {'acqTime'};
    kinStruct = getKinStruct(sucData,condFields,kinFields);
    
%% Get Condition Lists
    kinData = [Data.kinData];
    targetList = unique([trajStruct.target]);
    
%% Get Target Stats
    targetStatsStruct = struct('target',[],'total',[],'suc',[],'sucPct',[],'medianAcqTime',[]);
    %Get raw data
    targetData = [Data.targetData];
    targetID = [targetData.targetID];
    trialStatus = double([Data.trialStatus]);

    %Fill Struct
    structInd = 1;
    for target = targetList
        targetStatsStruct(structInd).target = target;
        %Total
        sucData = trialStatus([targetID==target]);
        totalTrials = size(sucData,2);
        sucTrials = sum(sucData);
        sucPct = 100*sucTrials/totalTrials;
        targetStatsStruct(structInd).total = totalTrials;
        targetStatsStruct(structInd).suc = sucTrials;
        targetStatsStruct(structInd).sucPct = sucPct;
        %Acq Time
        medianAcqTime = median([kinData([targetID==target]).acqTime]);
        targetStatsStruct(structInd).medianAcqTime = medianAcqTime;
        %Update Struct Ind
        structInd = structInd +1;
    end
    
%% Define color maps
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    cmap = customRainbow;
    hsvColor = rgb2hsv(cmap);
    hsvColor(:,2)=.3;
    cmap2 =hsv2rgb(hsvColor);
    
%% Plot results
    %Reach Trajectories and success Rates
    f = figure; f.Position = [10 10 1400 800];
    subplot(4,2,[1 3 5 7])
    for target = targetList
        tempData = trajStruct([trajStruct.target]==target);
        if ~isempty(tempData)
            numTrials = size(tempData.allDecoderCursorTraj,2);
            for trial = 1:round(numTrials/5):numTrials
                traj = tempData.allDecoderCursorTraj(trial).traj;
                plot(traj(:,1),traj(:,2),'Color',cmap2(target,:),'LineWidth',0.5)
                hold on
            end
            avgTraj = tempData.avgDecoderCursorTraj.traj;
            plot(avgTraj(:,1),avgTraj(:,2),'Color',cmap(target,:),'LineWidth',3);
        end
        targetAngle = (target-1)*45;
        %Get success info
        targetStats = targetStatsStruct([targetStatsStruct.target]==target);
        %Write stats on plot
        textDist = 100;
        textX = textDist*cosd(targetAngle);
        textY = textDist*sind(targetAngle);
            text(textX,textY,[num2str(round(targetStats.sucPct)),'% ','(',num2str(targetStats.suc),'/',num2str(targetStats.total),')'...
                newline,'MAT = ',num2str(round(targetStats.medianAcqTime))...
                ],'FontSize',10)
    end
    axis equal
    xlabel('x (mm)')
    ylabel('y (mm)')
    xlim([-100 100])
    ylim([-100 100])
    xticks([-100 100])
    yticks([-100 100])
    
    %Write stats on plot
    targetAngle = (target-1)*45;
    targetStats = targetStatsStruct([targetStatsStruct.target]==target);
    textDist = 100;
    textX = textDist*cosd(targetAngle);
    textY = textDist*sind(targetAngle);
        text(textX,textY,[num2str(round(targetStats.sucPct)),'% ','(',num2str(targetStats.suc),'/',num2str(targetStats.total),')'...
            newline,'MAT = ',num2str(round(targetStats.medianAcqTime))...
            ],'FontSize',10)
        
    %Forces
    for target = targetList
        tempData = trajStruct([trajStruct.target]==target);
        if ~isempty(tempData)
            numTrials = size(tempData.allForce,2);
            for trial = 1:round(numTrials/10):numTrials
                traj = tempData.allForce(trial).traj;
                time = tempData.allForce(trial).timestamps;
                for dir = 1:4
                    subplot(4,2,2*dir)
                    plot(time,traj(:,dir),'Color',cmap2(target,:),'LineWidth',0.5)
                    hold on
                end
            end
            avgTraj = tempData.avgForce.traj;
            avgTime = tempData.avgForce.timestamps;
            for dir = 1:4
                subplot(4,2,2*dir)
                plot(avgTime,avgTraj(:,dir),'Color',cmap(target,:),'LineWidth',3);
            end
        end

    end
    for dir = 1:4
        subplot(4,2,2*dir)
        xlabel('time (ms)')
        switch dir
            case 1
                ylabel('F_X (N)')
            case 2
                ylabel('F_Y (N)')
            case 3
                ylabel('F_Z (N)')
            case 4
                ylabel('F_{Total} (N)')
        end
    end
    sgtitle(dateStr)
    if saveFig
        saveas(gcf,[savePath,'reachTraj',dateStr,'.jpg'])
    end


% %% Title figure and Save
%     sgtitle(date)
%     if saveFig 
%         saveas(gcf,[savePath,'DCOdailyAnalysis',date,'.jpg'])
%     end