clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    task  = 'BCI';

    switch task
        case 'BCI'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\BCI Figures';
        case 'planning'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\Planning Figures';
        case 'iso'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\Iso Figures';
        case 'reaching'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\Reaching Figures';
        case 'multijoint BCI'
            saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Thesis Committee Meetings\Thesis Committee Update 1-16-22\Multijoint BCI Figures';
    end
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;
    
%% Load Data, subselect, get traj struct
    trajFields = {'smoothFR','marker'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;

    switch task
        case 'BCI'
        [Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
        [Data] = subselectForTrajDist(Data,task);
        trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
        condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
        %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-150},{'state','Step 2','first',0}};
        trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
        trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    end
    
%% Get posture and target lists
    trajStruct = trajStruct;
    switch task
        case {'BCI','planning','iso','multijoint BCI'}
            postureList = unique([trajStruct.posture]);
        case 'reaching'
            postureList = unique([fourPostureTrajStruct.posture]);
    end
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
    
%% For each posture, get hand positions
    resultStruct = struct('posture',[],'position',[]);
    structInd = 1;
    for posture = postureList
        tempData = trajStruct([trajStruct.posture]==posture);
        tempPos = [];
        for i = 1:size(tempData,2)
           tempPos = vertcat(tempPos,tempData(i).avgMarker.traj); 
        end
        resultStruct(structInd).posture = posture;
        resultStruct(structInd).position = mean(tempPos,1);
        structInd = structInd + 1;
    end

%% Get max dist
dist = vecnorm(resultStruct(5).position-resultStruct(1).position);