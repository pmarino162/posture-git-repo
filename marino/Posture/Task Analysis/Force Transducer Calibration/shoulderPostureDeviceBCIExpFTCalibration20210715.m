% clear; clc; clf; close all;

% %% Get bias in posture bci setup    
%     % Load Data
%     load('D:\Animals\Earl\2021\07\forceTransducerTesting20210714\blankBaseline\Earl20210714_blankBaseline_translated.mat')
% 
%     % Identify force transducer analog channels
%     algCh = Data(1).Definitions.analogChannelNames;
%     chName = {'Force X','Force Y','Force Z','Torque X','Torque Y','Torque Z'};
%     chIdx = NaN(size(chName));
%     for i = 1:length(chName)
%         if ~sum(strcmpi(algCh,chName{i}))
%             error('Force transducer channel %s not found',chName{i})
%         end
%         chIdx(i) = find(strcmpi(algCh,chName{i}),1);
%     end
%     
%     numTrials = size(Data,2);
%     biasData = [];
%     for trial = 1:numTrials
%        ftAlg = Data(trial).TrialData.analogData(:,chIdx);
%        biasData = vertcat(biasData,ftAlg);
%     end
%     
%     bias = mean(biasData,1);
%     
%% Examine force readouts with bias removed
%      load('D:\Animals\Earl\2021\07\forceTransducerTesting20210714\blankForces(foward-left-down)\Earl20210714_blankForces(foward-left-down)_translated.mat');
%     Data = getDataStruct(Data,'getForce',true,'forceSetup','shoulder_posture_bci');
    
%     load('\\actinium.univ.pitt.edu\batista\Animals\Earl\2021\07\forceTransducerTesting20210715\iso_force_config_autoMonkeyTest_LRFBUD\Earl20210715_iso_force_config_autoMonkeyTest_LRFBUD_translated.mat');
%     load('\\actinium.univ.pitt.edu\batista\Animals\Earl\2021\07\forceTransducerTesting20210715\iso_force_config_autoMonkeyTest_LRFBUD_level\Earl20210715_iso_force_config_autoMonkeyTest_LRFBUD_level_translated.mat');
%     load('\\actinium.univ.pitt.edu\batista\Animals\Earl\2021\08\elbowFTTest_LRUDFB\Earl20210809_elbowFTTest_LRUDFB_translated.mat');
    Data = getDataStruct(Data,'getForce',true,'forceSetup','shoulder_bar');
    allForce = [];
    allForceCursor = [];
    allAbsoluteForcesAndTorques = [];
    allTime = [];
    
    numTrials = size(Data,2);
    for trial = 1:numTrials
       allForce = vertcat(allForce,Data(trial).force.force);
       allForceCursor = vertcat(allForceCursor,Data(trial).force.forceCursor);
       allAbsoluteForcesAndTorques = vertcat(allAbsoluteForcesAndTorques,Data(trial).force.absoluteForcesAndTorques);
       time = Data(trial).force.time;
       if trial > 1
           time = time + allTime(end);
       end
       allTime = [allTime,time];
    end
    
%% Plot results
    %Force 
    f = figure; f.Position = [10 10 1400 800];   
    for dir = 1:4
        subplot(4,1,dir)
        plot(allTime,allForce(:,dir));
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
    sgtitle('MonkeyHost Forces')
    
    %Force Cursor
    f = figure; f.Position = [10 10 1400 800];   
    for dir = 1:3
        subplot(3,1,dir)
        plot(allTime,allForceCursor(:,dir));
        xlabel('time (ms)')
        switch dir
            case 1
                ylabel('x (mm)')
            case 2
                ylabel('y (mm)')
            case 3
                ylabel('z (mm)')
        end
    end
    sgtitle('Force Cursor')
    
    %Absolute Forces and Torques
    f = figure; f.Position = [10 10 1400 800];   
    for dir = 1:6
        subplot(6,1,dir)
        plot(allTime,allAbsoluteForcesAndTorques(:,dir));
        xlabel('time (ms)')
        switch dir
            case 1
                ylabel('F_X (N)')
            case 2
                ylabel('F_Y (N)')
            case 3
                ylabel('F_Z (N)')
            case 4
                ylabel('T_X (N.m)')
            case 5
                ylabel('T_Y (N.m)')
            case 6
                ylabel('T_Z (N.m)') 
        end
    end
    sgtitle('Absolute Forces And Torques')