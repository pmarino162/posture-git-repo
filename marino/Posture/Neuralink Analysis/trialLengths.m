clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\Figure 3\Projections';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Run loop for each dataset      
    %Load data
    dataset = 'E20200316';
    [Data,zScoreParams] = loadData(dataset);

%% Plot move time
    moveTime = [Data.kinData];
    moveTime = [moveTime.moveTime];
    figure;
    histogram(moveTime);
    xlabel('Move Time (ms)')
    ylabel('Count')

