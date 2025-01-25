clear; clc; clf; close all;

%% Load, Preprocess, and Label Data    
    %8/21/2019
    date = '20190821';
    [Data,GPFAParams,DecoderParams] = loadEarlData20190821;
 
%% Create Task ID 
    taskIDs = struct('ID',[],'Task',[]);
    taskIDs(1).ID = 1; taskIDs(1).Task = 'BC';
    taskIDs(2).ID = 2; taskIDs(2).Task = 'HC';
    taskIDs(3).ID = 3; taskIDs(3).Task = 'Iso';

%% Take Decode Plane Condition Averages 
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    cmap = customRainbow;

    BC = Data(strcmpi({Data.Task},'BC'));
        condFields = {{'target','targetData','targetID'}};
        trajFields = {'decoderWorkTraj'};
        trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'CenterOutCenter_BC_TouchBar'};
        trialInclStates(1).inclStates = {'Step 1'}; trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {-50};
        trialInclStates(1).addTimeToEnd = {0};
        trajStruct = getTrajStruct(BC,condFields,trajFields,trialInclStates);
        
        
    figure
    for target = 1:8
        traj = trajStruct([trajStruct.target]==target).avgDecoderWorkTraj.traj;
        plot(traj(:,1),traj(:,2))
        hold on
    end
    
    
    xlabel('projection onto decoder X (mm)')
    ylabel('projection on decoder Y (mm)')
    title('Brain Control')
    
    HC = Data(strcmpi({Data.Task},'HC'));
        condFields = {{'target','targetData','targetID'}};
        trajFields = {'decoderWorkTraj','marker'};
        trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'HC_CenterOut'};
        trialInclStates(1).inclStates = {'Step 2'}; trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {-50};
        trialInclStates(1).addTimeToEnd = {0};
        trajStruct = getTrajStruct(HC,condFields,trajFields,trialInclStates);
        
        
    figure
    for target = 1:8
        traj = trajStruct([trajStruct.target]==target).avgDecoderWorkTraj.traj;
        plot(traj(:,1),traj(:,2))
        hold on
    end
    xlabel('projection onto decoder X (mm)')
    ylabel('projection on decoder Y (mm)')
    title('Reaching')

    
    Iso = Data(strcmpi({Data.Task},'Iso'));
        condFields = {{'target','targetData','targetID'}};
        trajFields = {'decoderWorkTraj'};
        trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'IsometricForce_1D'};
        trialInclStates(1).inclStates = {'Target'}; trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {-50};
        trialInclStates(1).addTimeToEnd = {0};
        trajStruct = getTrajStruct(Iso,condFields,trajFields,trialInclStates);
        
        
    figure
    for target = [3,7]
        traj = trajStruct([trajStruct.target]==target).avgDecoderWorkTraj.traj;
        plot(traj(:,1),traj(:,2))
        hold on
    end
    
    xlabel('projection onto decoder X (mm)')
    ylabel('projection on decoder Y (mm)')
    title('Isometric Force')
    
%% Plot
% % Load Colormaps 
%     load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
%     cmap = customRainbow;
% 
% 
% figure
% for target = 1:8
%     traj = trajStruct([trajStruct.target]==target).avgDecoderWorkTraj.traj;
%     plot(traj(:,1),traj(:,2))
%     hold on
% end