
%% BCI task

load('D:\Animals\Earl\2020\03\20200316\05_E30_grid\Earl20200316_05_E30_grid_SI_translated.mat')

forceWindow = Data(1).Parameters.StateTable(4).Hand.window(2,:);
forceBounds = forceWindow(5);
 [Data] = getForceAndForceCursor(Data,setup)
bciTargetLocations

%% Isometric force task
isoTargetLocations