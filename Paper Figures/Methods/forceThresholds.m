clear; clc

%% Monkey E BCI calibration 
    load('D:\Animals\Earl\2020\03\20200316\04_E30_gradualTraining\Earl20200316_04_E30_gradualTraining_SI_translated.mat');
    taskForceCursorWindow = Data(1).Parameters.StateTable(5).Hand.window(2,:);
    forceCursorRange = taskForceCursorWindow(5);
    scaleFactor = Data(1).Definitions.forceTransformation.scaling(2);
    maxForce = (forceCursorRange/2)/scaleFactor;

%% Monkey E BCI task
    load('D:\Animals\Earl\2020\03\20200316\05_E30_grid\Earl20200316_05_E30_grid_SI_translated.mat')
    taskForceCursorWindow = Data(1).Parameters.StateTable(5).Hand.window(2,:);
    forceCursorRange = taskForceCursorWindow(5);
    scaleFactor = Data(1).Definitions.forceTransformation.scaling(2);
    maxForce = (forceCursorRange/2)/scaleFactor;

%% Monkey E Isometric force task
    load('D:\Animals\Earl\2020\01\20200116\E20200116.mat')
    targetCursorWindow = Data(1).targetData.targetSize(2);
    forceCursorRange = Data(2).targetData.targetLoc(2)-Data(1).targetData.targetLoc(2);
    load('D:\Animals\Earl\2020\01\20200116\03_isoForceE15\Earl20200116_03_isoForceE15_SI_SW_translated.mat')
    scaleFactor = Data(1).Definitions.forceTransformation.scaling(2);
    maxForce = (forceCursorRange/2)/scaleFactor;
    targetForceWindow = (targetCursorWindow/2)/scaleFactor;
    
%% Monkey E Multiple tasks BCI
    load('D:\Animals\Earl\2020\03\20200314\04_neutral_grid\Earl20200314_04_neutral_grid_SI_translated.mat')
    taskForceCursorWindow = Data(1).Parameters.StateTable(5).Hand.window(2,:);
    forceCursorRange = taskForceCursorWindow(5);    
    scaleFactor = Data(1).Definitions.forceTransformation.scaling(2);
    maxForce = (forceCursorRange/2)/scaleFactor;

%% Monkey E Multiple Tasks Iso
    load('D:\Animals\Earl\2020\03\20200314\E20200314.mat')   
    targetCursorWindow = Data(294).targetData.targetSize(2);
    forceCursorRange = Data(294).targetData.targetLoc(2)-Data(293).targetData.targetLoc(2);
    load('D:\Animals\Earl\2020\03\20200314\06_neutral_isoForce\Earl20200314_06_neutral_isoForce_SI_translated.mat')
    scaleFactor = Data(1).Definitions.forceTransformation.scaling(2);
    maxForce = (forceCursorRange/2)/scaleFactor;
    targetForceWindow = targetCursorWindow/scaleFactor;
