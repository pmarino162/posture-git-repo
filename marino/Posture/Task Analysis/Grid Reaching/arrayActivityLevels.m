clear; clc; clf; close all;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = rainbow2d;
    
%% Load Data 
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
    exclCh = [];%[44 87 88 77 78 71 67 69 118]
    Data = getDataStruct20211210(Data,'getSpikes',true,'getSorts',false,'exclCh',[],'exclZero',false,'getMarker',true,'centerMarker',false,'getKin',true);
    [Data] = cleanData20210706(Data);
    Data = Data([Data.trialStatus]==1);
    [Data,postureIDs] = labelPostures20210706(Data);
%     allTrialPostures = [Data.conditionData]; allTrialPostures = [allTrialPostures.postureID];
%     Data = Data(ismember(allTrialPostures,[1,2,4,5]));
    

%% Get trajStruct 
    trajFields = {'smoothFR'};
    binWidth = 1;
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'status','trialStatus'}};
    trialInclStates(1).trialName = {'GridReaching'};
    %Reach traj struct
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-100},{'state','Target Hold','first',0}};    
    reachTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    %Delay traj Struct
    trialInclStates(1).inclStates = {{'state','Delay','first',-100},{'state','Target Acquire','first',0}};
    delayTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%% z-score reach neural activity, plot 
    %Preallocate 
    numTraj = size(reachTrajStruct.allSmoothFR,2);
    numObs = 0;
    numDims = size(reachTrajStruct.allSmoothFR(1).traj,2);
    for i = 1:numTraj
        numObs = length(reachTrajStruct.allSmoothFR(i).timestamps)+ numObs;
    end
    allNeuralData = NaN(numObs,numDims);
    j = 1;
    for i = 1:numTraj
        traj = reachTrajStruct.allSmoothFR(i).traj;
        numSteps = size(traj,1);
        allNeuralData(j:j+numSteps-1,:) = traj;
        j = j+numSteps;
    end
    chMeans = mean(allNeuralData,1);
    chStd = std(allNeuralData);
    
% Plot colormap
    FRtraj = reachTrajStruct.avgSmoothFR.traj;
    FRtrajZ = (FRtraj-chMeans)./chStd;
    
    figure
    imagesc(FRtraj')
    xlabel('time (ms)')
    ylabel('channel')

    figure
    imagesc(FRtrajZ')
    xlabel('time (ms)')
    ylabel('channel')

    figure
    bar(chMeans)
    xlabel('channel')
    ylabel('mean FR (Hz)')
    
    
%% z-score delay neural activity, plot 
    %Preallocate 
    numTraj = size(delayTrajStruct.allSmoothFR,2);
    numObs = 0;
    numDims = size(delayTrajStruct.allSmoothFR(1).traj,2);
    for i = 1:numTraj
        numObs = length(delayTrajStruct.allSmoothFR(i).timestamps)+ numObs;
    end
    allNeuralData = NaN(numObs,numDims);
    j = 1;
    for i = 1:numTraj
        traj = delayTrajStruct.allSmoothFR(i).traj;
        numSteps = size(traj,1);
        allNeuralData(j:j+numSteps-1,:) = traj;
        j = j+numSteps;
    end
    chMeans = mean(allNeuralData,1);
    chStd = std(allNeuralData);
    
% Plot colormap
    FRtraj = delayTrajStruct.avgSmoothFR.traj;
    FRtrajZ = (FRtraj-chMeans)./chStd;
    
    figure
    imagesc(FRtraj')
    xlabel('time (ms)')
    ylabel('channel')

    figure
    imagesc(FRtrajZ')
    xlabel('time (ms)')
    ylabel('channel')

    figure
    bar(chMeans)
    xlabel('channel')
    ylabel('mean FR (Hz)')
    