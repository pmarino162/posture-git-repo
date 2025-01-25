clear
clc
clf
close all

%% Load Data
Data = load('D:\Animals\Earl\2021\05\20210517\Internal\Earl20210517_Internal_translated.mat');
%Change Z scale 
    for trial = 1:size(Data.Data,2)
       Data.Data(trial).Definitions.forceTransformation.scaling(1,3) = 25;
    end
Data = getDataStruct(Data.Data,'getForce',true,'forceSetup','shoulder_bar'); 
Data = Data([Data.trialStatus]==1);

%% Create trajStruct
    condFields = {{'target','targetData','targetID'}};
    trajFields = {'force','forceCursor'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'IsometricForce_1D'};
        trialInclStates(1).inclStates = {'Target','Target Hold'}; trialInclStates(1).inclOccurrence = {'last','last'};
        trialInclStates(1).addTimeToBeginning = {-100,0};
        trialInclStates(1).addTimeToEnd = {0,0};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
    
%% Plot Forces 
for target = [1,3,5,7]
    %Forces 
   tempData = trajStruct([trajStruct.target]==target).avgForce;
   figure
   hold on
   for plotInd = 1:4
        subplot(4,1,plotInd)
        force = tempData.traj(:,plotInd);
        time = tempData.timestamps;
        plot(time,force)
        xlabel('time (ms)')
        switch plotInd
            case 1
                ylabel('Fx (N)')
            case 2
                ylabel('Fy (N)')
            case 3
                ylabel('Fz (N)')
            case 4
                ylabel('Total Force (N)')
        end
   end
   sgtitle(['Forces: Target ',num2str(target)])
   
   %Force Cursor
   tempData = trajStruct([trajStruct.target]==target).avgForceCursor;
   figure
   hold on
   for plotInd = 1:3
        subplot(3,1,plotInd)
        forceCursor = tempData.traj(:,plotInd);
        time = tempData.timestamps;
        plot(time,forceCursor)
        xlabel('time (ms)')
        switch plotInd
            case 1
                ylabel('x (mm)')
            case 2
                ylabel('y (mm)')
            case 3
                ylabel('z (mm)')
        end
   end
   sgtitle(['Force Cursor: Target ',num2str(target)])
   
end