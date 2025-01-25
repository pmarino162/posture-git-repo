clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20221123\Figures\Figure 2\Alignment Indices\20230125\fullSpace';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = orli;
    scmap = customRainbow;
    
%% Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],'modelR2',[],'modelError',[],'noShiftR2',[],'noShiftError',[],'shuffleR2',[],'shuffleError',[],'selfR2',[],'selfError',[]);
    structInd = 1;
    
%% Get trajStruct
    dataset = 'E20200316';
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
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

%% Get timestamps, posture & target lists
 numTimestamps = [];
    numCondTraj = [];
    for i = 1:size(trajStruct,2)
       numTraj = size(trajStruct(i).allSmoothFR,2);
       numCondTraj = [numCondTraj,numTraj];
       for j = 1:numTraj
          numTimestamps = [numTimestamps,size(trajStruct(i).allSmoothFR(j).timestamps,2)]; 
       end
    end
    figure
    histogram(numTimestamps)
    xlabel('Number of 25ms bins')
    ylabel('Number of trials')
    
    figure
    histogram(numCondTraj)
    xlabel('Number of trials')
    ylabel('Number of conditions')
    
    % Get minimum number of trials and timestamps
    [minNumTimestamps,i] = min(numTimestamps);
    [minNumCondTraj,i] = min(numCondTraj);
    
    %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);
    
%% Do computation
    %Parameters
    numIterations = 100;
    numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw
    numPts = minNumTimestamps; %Number of points from each trajectory

    for i = 1:numIterations 
        %Create 2 sets of condition averages
        trajStruct1 = struct('posture',[],'target',[],'traj',[]);
        trajStruct2 = struct('posture',[],'target',[],'traj',[]);
        for j = 1:size(trajStruct,2)
            numTraj = size(trajStruct(j).allSmoothFR,2);
            sampInd1 = randsample(numTraj,numSample);
            sampInd2 = 
            traj1 = getAvgTraj20211210(trajStruct(j).allSmoothFR(sampInd1),binWidth);
            traj2 = getAvgTraj20211210(trajStruct(j).allSmoothFR(sampInd2),binWidth);
            
            posture = trajStruct(j).posture; target = trajStruct(j).target;
            
            trajStruct1(j).posture = posture; 
            trajStruct1(j).target = target;
            trajStruct1(j).traj = traj1;
            trajStruct2(j).posture = posture; 
            trajStruct2(j).target = target;
            trajStruct2(j).traj = traj2;
        end
        
        %Get model, no shift, and shuffle predictions
        modelStruct = struct('posture',[],'target',[],'traj',[]);
        noShiftStruct = struct('posture',[],'target',[],'traj',[]);
        shuffleStruct = struct('posture',[],'target',[],'traj',[]);
        for posture = postureList
           for target = targetList 
                %Get target shape prediction - average over posture for
                %target (excluding held out traj)
                
                %Get posture shift prediction - repeat procedure for other
                %targets, and figure out how to shift those to get to held
                %out posture
           end
        end
        
        %For each condition, make each type of prediction; assess error
        
    end
    
    

%% Plot results
     
     