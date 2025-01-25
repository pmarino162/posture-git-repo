clear; clc; clf; close all; 

%% Load Data 
    %Earl - 3/17/2020
    date = '20200317';
    exclCh =  [3 12 16 20 31 73 85 92 94];
    binWidth = 1;
    [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
     N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200317(binWidth);
    numPostures = 5;
    numTargets = 8;
    
%% For each trial, get internal, neutral, and external projection
    numTrials = size(Data,2);
    for trial = 1:numTrials
        rawSpikeBins = Data(trial).Decoder.rawSpikeBins;
        %N00
        N00Decoder = struct('rawSpikeBins',rawSpikeBins,'Parameters',N00DecoderParams);
        [N00GPFATraj,N00WorkTraj,N00NullTraj] = getDecoderGPFATraj(N00Decoder);
        N00WorkTraj(:,2) = -N00WorkTraj(:,2);
        Data(trial).Decoder.N00WorkTraj = N00WorkTraj;
        %I30
        I30Decoder = struct('rawSpikeBins',rawSpikeBins,'Parameters',I30DecoderParams);
        [I30GPFATraj,I30WorkTraj,I30NullTraj] = getDecoderGPFATraj(I30Decoder);
        I30WorkTraj(:,2) = -I30WorkTraj(:,2);
        Data(trial).Decoder.I30WorkTraj = I30WorkTraj;
        %E30
        E30Decoder = struct('rawSpikeBins',rawSpikeBins,'Parameters',E30DecoderParams);
        [E30GPFATraj,E30WorkTraj,E30NullTraj] = getDecoderGPFATraj(E30Decoder);
        E30WorkTraj(:,2) = -E30WorkTraj(:,2);
        Data(trial).Decoder.E30WorkTraj = E30WorkTraj;
    end

%% Create trajStruct
    condFields = {{'posture','postureData','postureID'},{'target','targetData','target1ID'}};
    trajFields = {'N00WorkTraj','I30WorkTraj','E30WorkTraj'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
        trialInclStates(1).inclStates = {'Step 1'}; trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {-100};
        trialInclStates(1).addTimeToEnd = {0};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);

%% Plot Trajectories
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    %N00 Decoder    
    figure 
    xDim = 1; yDim = 2;
    for posture = [1,3,5]
        for target = 1:numTargets
            traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgN00WorkTraj.traj;
            plot(traj(:,xDim),traj(:,yDim),'Color',cmap(posture,:),'LineWidth',2); 
            hold on
        end
    end
    title('N00 Decoder')
    
%E30 Decoder
figure 
xDim = 1; yDim = 2;
for posture = [1,5]
    for target = 1:numTargets
        traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgE30WorkTraj.traj;
        plot(traj(:,xDim),traj(:,yDim),'Color',cmap(posture,:),'LineWidth',2); 
        hold on
    end
end
title('E30 Decoder')

%I30 Decoder
figure 
xDim = 1; yDim = 2;
for posture = [1,5]
    for target = 1:numTargets
        traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).avgI30WorkTraj.traj;
        plot(traj(:,xDim),traj(:,yDim),'Color',cmap(posture,:),'LineWidth',2); 
        hold on
    end
end
title('I30 Decoder')